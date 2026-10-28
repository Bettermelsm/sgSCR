#!/usr/bin/env bash
# run_step.sh — Step 调度与环境切换（主环境 Rscript / 隔离环境 conda run）
# Validates: Requirements 4.1, 4.2, 4.4, 4.5（design.md「C. 环境隔离架构」）
#
# 用法（远程执行，在 remote_work_dir 下）:
#   bash scripts/common/run_step.sh <step_id> --config config/config.yaml
#   bash scripts/common/run_step.sh --from s3 --config config/config.yaml   # 断点续跑
#
# 环境查表见 scripts/common/env_resolver.R 的 .sgs_STEP_ENV（单一权威来源，
# shell 侧重复此表，两处须同步维护——新增 Step 时两处都要登记）。

set -euo pipefail
die() { echo "[ERROR] [run_step] $*" >&2; exit 1; }

declare -A STEP_ENV=(
  [s0]="main" [s1]="main" [s2]="main" [s2b]="main" [s3]="main" [s4]="main"
  [s5]="main" [s6]="main" [s7]="cellchat" [s8]="trajectory" [s9]="scenic"
  [s10]="main" [s10b]="main" [s11]="main"
)

STEP=""; FROM=""; CONFIG=""
while [ $# -gt 0 ]; do
  case "$1" in
    --config) CONFIG="$2"; shift 2 ;;
    --from)   FROM="$2"; shift 2 ;;
    --)       shift; break ;;
    -*)       die "unknown flag: $1" ;;
    *)        STEP="$1"; shift ;;
  esac
done
[ -n "$CONFIG" ] && [ -f "$CONFIG" ] || die "--config <file> required"
[ -n "$STEP" ] || [ -n "$FROM" ] || die "need <step_id> or --from <step_id>"
[ -z "$STEP" ] && STEP="$FROM"

# P12-a A2：remote.mode=local 时本脚本本机直接执行（无需远程封装）；
# =ssh 时本脚本预期已在 remote_work_dir 内被调用，此处仅提示而非拦截
REMOTE_MODE="$(python3 - "$CONFIG" <<'PY'
import sys, yaml
cfg = yaml.safe_load(open(sys.argv[1], encoding="utf-8")) or {}
print(((cfg.get("remote") or {}).get("mode")) or "local")
PY
)"
case "$REMOTE_MODE" in
  local) : ;;
  ssh)   echo "[INFO] [run_step] remote.mode=ssh — 确认本脚本运行在 remote_work_dir 内（remote_exec.sh sync/exec 启动）" >&2 ;;
  *)     die "remote.mode 必须为 local 或 ssh，当前: '$REMOTE_MODE'" ;;
esac

# P12-a D6：config.run_id 与输出目录 .run_id 标记一致性校验
# （run1 误写事故固化：目录错配 → 立即拒绝，绝不静默写入别的 run）
CFG_RUN_ID="$(python3 - "$CONFIG" <<'PY'
import sys, yaml
cfg = yaml.safe_load(open(sys.argv[1], encoding="utf-8")) or {}
print((cfg.get("run") or {}).get("run_id") or "")
PY
)"
MANIFEST_PATH="$(python3 - "$CONFIG" <<'PY'
import sys, yaml
cfg = yaml.safe_load(open(sys.argv[1], encoding="utf-8")) or {}
print(((cfg.get("paths") or {}).get("manifest_path")) or "")
PY
)"
if [ -n "$MANIFEST_PATH" ] && [ -f "$(dirname "$MANIFEST_PATH")/.run_id" ]; then
  DIR_RUN_ID="$(cat "$(dirname "$MANIFEST_PATH")/.run_id")"
  [ "$DIR_RUN_ID" = "$CFG_RUN_ID" ] \
    || die "run_id 不一致: config.run_id='$CFG_RUN_ID' 但输出目录 .run_id='$DIR_RUN_ID' — 请核对 --config 指向（run1 误写事故防护）"
fi

KEY="${STEP%%_*}"   # "s7_cellchat" -> "s7"
[ -n "${STEP_ENV[$KEY]:-}" ] || die "unregistered step: $STEP (known keys: ${!STEP_ENV[*]})"
ENV_KEY="${STEP_ENV[$KEY]}"

# 从 config 读取该 step 对应的 conda env 名
ENV_NAME="$(python3 - "$CONFIG" "$ENV_KEY" <<'PY'
import sys, yaml
cfg = yaml.safe_load(open(sys.argv[1], encoding="utf-8"))
print((cfg.get("envs") or {}).get(sys.argv[2]) or "")
PY
)"
[ -n "$ENV_NAME" ] || die "envs.$ENV_KEY not set in $CONFIG"

SCRIPT="$( { ls scripts/${KEY}_*.R scripts/${KEY}_*.py 2>/dev/null || true; } | head -1)"
[ -n "$SCRIPT" ] || die "no script matched scripts/${KEY}_*.R|.py"
case "$SCRIPT" in
  *.py) INTERP="python" ;;
  *)    INTERP="Rscript" ;;
esac

RESOLVED=()
[ -n "$FROM" ] && RESOLVED=(--from "$FROM")

# 启动前校验环境存在；main 允许两种形态：系统 R（不在 conda 里→直接 Rscript）
# 或 conda env（如 sgscr_main → conda run）。隔离环境必须是 conda env。
CONDA_BIN="$(command -v conda || echo "$HOME/anaconda3/bin/conda")"
env_in_conda() {
  [ -x "$CONDA_BIN" ] && "$CONDA_BIN" env list 2>/dev/null | grep -qE "(^|/| )${1}\s*$"
}

if [ "$ENV_KEY" = "main" ] && ! env_in_conda "$ENV_NAME"; then
  echo "[INFO] [run_step] $STEP -> system $INTERP: $SCRIPT ${RESOLVED[*]:-}"
  exec "$INTERP" "$SCRIPT" --config "$CONFIG" "${RESOLVED[@]}"
else
  env_in_conda "$ENV_NAME" \
    || die "conda env '$ENV_NAME' (envs.$ENV_KEY) not found — create it first (scripts/common/setup_envs.sh)"
  echo "[INFO] [run_step] $STEP -> conda env $ENV_NAME ($INTERP): $SCRIPT ${RESOLVED[*]:-}"
  exec "$CONDA_BIN" run -n "$ENV_NAME" "$INTERP" "$SCRIPT" --config "$CONFIG" "${RESOLVED[@]}"
fi
