#!/usr/bin/env bash
# remote_exec.sh — 远程连接容错 + 脚本/配置同步（本机侧运行）
# Validates: Requirements 3.1, 3.2 (design.md「D. 远程执行机制」)
#
# 用法:
#   remote_exec.sh sync  --config <config.yaml>                # 同步 scripts/ 与 config/ 到远程
#   remote_exec.sh exec  --config <config.yaml> -- <cmd...>    # 在远程执行命令
#   remote_exec.sh fetch --config <config.yaml> -- <remote> <local>  # 取回产物
#
# 连接策略：优先 ssh 直连（remote.ssh_target，可为 ssh config Host 别名）；
# 失败时尝试 frp 隧道（remote.frp_fallback.host + port）；
# 两者皆不可用则终止并分别输出两种方式的失败详情。

set -uo pipefail
SKILL_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"

log()  { echo "[INFO] [remote_exec] $*" >&2; }
warn() { echo "[WARN] [remote_exec] $*" >&2; }
die()  { echo "[ERROR] [remote_exec] $*" >&2; exit 1; }

PY="$(command -v python3 || command -v python)"
[ -n "$PY" ] || die "python not found on PATH (needed for lightweight YAML read)"

# 极简 YAML 标量读取（仅 shell 侧引导；R 侧用 yaml 包完整解析）
yaml_get() { # yaml_get <file> <dot.path>
  "$PY" - "$1" "$2" <<'PY'
import sys, yaml
try:
    cfg = yaml.safe_load(open(sys.argv[1], encoding="utf-8"))
except Exception:
    print(""); sys.exit(0)
cur = cfg
for p in sys.argv[2].split("."):
    cur = cur.get(p) if isinstance(cur, dict) else None
    if cur is None: break
print(cur if cur is not None else "")
PY
}

CONFIG=""; MODE="${1:-}"; shift || true
while [ $# -gt 0 ]; do
  case "$1" in
    --config) CONFIG="$2"; shift 2 ;;
    --) shift; break ;;
    *) shift ;;
  esac
done
[ -n "$CONFIG" ] && [ -f "$CONFIG" ] || die "--config <existing file> required"

SSH_TARGET="$(yaml_get "$CONFIG" remote.ssh_target)"
# P12-a A2：remote.mode=local → 单机模式，sync 无需同步、exec/fetch 本机直做
REMOTE_MODE="$(yaml_get "$CONFIG" remote.mode)"
case "$REMOTE_MODE" in
  ""|local|ssh) ;;
  *) die "remote.mode 必须为 local 或 ssh，当前: '$REMOTE_MODE'" ;;
esac
if [ "$REMOTE_MODE" = "local" ]; then
  case "$MODE" in
    sync)  log "remote.mode=local：单机模式，无需同步"; exit 0 ;;
    exec)  [ $# -gt 0 ] || die "exec requires -- <command...>"; log "local 模式：本机执行: $*"; exec "$@" ;;
    fetch) REMOTE_P="${1:-}"; LOCAL_P="${2:-}"
           [ -n "$REMOTE_P" ] && [ -n "$LOCAL_P" ] || die "fetch requires -- <remote_path> <local_path>"
           cp -r "$REMOTE_P" "$LOCAL_P" || die "fetch failed: $REMOTE_P"
           log "local 模式：copied $REMOTE_P -> $LOCAL_P"; exit 0 ;;
    *)     die "unknown mode: '$MODE' (sync|exec|fetch)" ;;
  esac
fi
FRP_HOST="$(yaml_get "$CONFIG" remote.frp_fallback.host)"
FRP_PORT="$(yaml_get "$CONFIG" remote.frp_fallback.port)"
WORK_DIR="$(yaml_get "$CONFIG" remote.remote_work_dir)"
[ -n "$SSH_TARGET" ] || die "remote.ssh_target not set in config"

# 容错远程执行：直连优先，frp 兜底；两者失败详情都输出
remote_run() {
  local out
  if out=$(ssh -o ConnectTimeout=8 -o BatchMode=yes "$SSH_TARGET" "$*" 2>&1); then
    echo "$out"; return 0
  fi
  warn "ssh direct to '$SSH_TARGET' failed:"
  echo "$out" >&2
  if [ -n "$FRP_HOST" ] && [ -n "$FRP_PORT" ]; then
    if out=$(ssh -o ConnectTimeout=8 -o BatchMode=yes -p "$FRP_PORT" "$FRP_HOST" "$*" 2>&1); then
      log "succeeded via frp tunnel $FRP_HOST:$FRP_PORT"
      echo "$out"; return 0
    fi
    warn "frp tunnel $FRP_HOST:$FRP_PORT also failed:"
    echo "$out" >&2
  else
    warn "frp_fallback not configured in config"
  fi
  die "both ssh direct and frp tunnel unavailable (details above)"
}

case "$MODE" in
  sync)
    [ -n "$WORK_DIR" ] || die "remote.remote_work_dir not set in config"
    remote_run "mkdir -p '$WORK_DIR'"
    # 先试直连 scp，再试 frp scp（rsync 远端若无则 scp 兜底，保持简单）
    if scp -q -r "$SKILL_DIR/scripts" "$SKILL_DIR/config" "$SSH_TARGET:$WORK_DIR/"; then
      log "synced scripts/ + config/ -> $SSH_TARGET:$WORK_DIR (direct)"
    elif [ -n "$FRP_HOST" ] && [ -n "$FRP_PORT" ] && \
         scp -q -r -P "$FRP_PORT" "$SKILL_DIR/scripts" "$SKILL_DIR/config" "$FRP_HOST:$WORK_DIR/"; then
      log "synced scripts/ + config/ -> $FRP_HOST:$WORK_DIR (frp)"
    else
      die "sync failed via both direct scp and frp scp"
    fi
    ;;
  exec)
    [ $# -gt 0 ] || die "exec requires -- <command...>"
    remote_run "$@"
    ;;
  fetch)
    REMOTE_P="${1:-}"; LOCAL_P="${2:-}"
    [ -n "$REMOTE_P" ] && [ -n "$LOCAL_P" ] || die "fetch requires -- <remote_path> <local_path>"
    scp -q "$SSH_TARGET:$REMOTE_P" "$LOCAL_P" || die "fetch failed: $REMOTE_P"
    log "fetched $REMOTE_P -> $LOCAL_P"
    ;;
  *) die "unknown mode: '$MODE' (sync|exec|fetch)" ;;
esac
