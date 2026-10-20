#!/usr/bin/env bash
# setup_envs.sh — 创建三个隔离 conda 环境（SGlcl01 上运行，mamba 优先）
# Validates: Requirements 4.1
#
# 用法: bash scripts/common/setup_envs.sh [--main-pkgs-only]
# 环境名以 config 的 envs.* 为准；本脚本用设计默认名，如改过 config 请传参覆盖：
#   bash setup_envs.sh --cellchat <name> --trajectory <name> --scenic <name>
#
# 说明：
# - main 环境（Seurat 等）通常 SGlcl01 已有，默认跳过创建，只装缺失的轻量包
# - Monocle3 / R 版 SCENIC 是 Config 可选项，此处不装（迭代项再评估）
# - pySCENIC 走 conda-forge + pip（pyscenic 在 pip 上更新）

set -euo pipefail
die() { echo "[ERROR] [setup_envs] $*" >&2; exit 1; }
info() { echo "[INFO] [setup_envs] $*" >&2; }

# 非交互 shell PATH 常无 conda/mamba，按常见安装位置探测
find_installer() {
  local c
  for c in mamba conda "$HOME/miniconda3/bin/mamba" "$HOME/anaconda3/bin/mamba" \
           "$HOME/miniconda3/bin/conda" "$HOME/anaconda3/bin/conda" \
           /opt/miniconda3/bin/conda /opt/anaconda3/bin/conda; do
    if [ "$c" = "mamba" ] || [ "$c" = "conda" ]; then
      command -v "$c" >/dev/null 2>&1 && { echo "$c"; return 0; }
    elif [ -x "$c" ]; then
      echo "$c"; return 0
    fi
  done
  return 1
}
INSTALLER="$(find_installer)" || die "conda/mamba not found in PATH or common locations"

CC_ENV="sgscr_cellchat"; TJ_ENV="sgscr_trajectory"; SC_ENV="sgscr_scenic"
MAIN_ENV=""   # 传 --main <name> 时创建主环境（服务器无现成 Seurat 时使用）
while [ $# -gt 0 ]; do
  case "$1" in
    --cellchat)   CC_ENV="$2"; shift 2 ;;
    --trajectory) TJ_ENV="$2"; shift 2 ;;
    --scenic)     SC_ENV="$2"; shift 2 ;;
    --main)       MAIN_ENV="$2"; shift 2 ;;
    *) die "unknown arg: $1" ;;
  esac
done

env_exists() { "$INSTALLER" env list 2>/dev/null | grep -qE "(^|/)$1\s*$"; }

ensure_env() { env_exists "$1" && { info "env exists, skip: $1"; return 0; } || true; }

# ---- sgscr_main（可选）：R 主分析环境（服务器无现成 Seurat 时） ----
if [ -n "$MAIN_ENV" ]; then
  if ! env_exists "$MAIN_ENV"; then
    info "creating $MAIN_ENV (R + Seurat/SingleR/clustree)..."
    "$INSTALLER" create -y -n "$MAIN_ENV" -c conda-forge       "r-base>=4.3" r-seurat bioconductor-singler r-clustree       r-yaml r-jsonlite r-testthat r-hedgehog r-future.apply
  else
    info "$MAIN_ENV exists; verifying Seurat loadable"
    "$INSTALLER" run -n "$MAIN_ENV" Rscript -e 'stopifnot(requireNamespace("Seurat", quietly=TRUE))'       || die "Seurat not loadable in $MAIN_ENV"
  fi
fi

# ---- sgscr_cellchat: R + CellChat + NMF/circlize ----
# conda-forge 无 r-cellchat（2026-10-02 实测 solve 失败）：先建基础 env，再 BiocManager 装
if ! env_exists "$CC_ENV"; then
  info "creating $CC_ENV (R base + deps; CellChat via BiocManager)..."
  "$INSTALLER" create -y -n "$CC_ENV" -c conda-forge \
    "r-base>=4.3" r-seurat r-nmf r-circlize r-future.apply r-remotes \
    r-biocmanager
  "$INSTALLER" run -n "$CC_ENV" Rscript -e \
    'BiocManager::install("CellChat", ask=FALSE, update=FALSE)' \
    || die "BiocManager install CellChat failed in $CC_ENV"
fi
info "verifying CellChat loadable in $CC_ENV"
"$INSTALLER" run -n "$CC_ENV" Rscript -e 'stopifnot(requireNamespace("CellChat", quietly=TRUE))' \
  || die "CellChat not loadable in $CC_ENV"

# ---- sgscr_trajectory: R + slingshot ----
if ! env_exists "$TJ_ENV"; then
  info "creating $TJ_ENV (R + slingshot)..."
  "$INSTALLER" create -y -n "$TJ_ENV" -c conda-forge \
    "r-base>=4.3" bioconductor-slingshot r-seurat
else
  info "$TJ_ENV exists; verifying slingshot loadable"
  "$INSTALLER" run -n "$TJ_ENV" Rscript -e 'stopifnot(requireNamespace("slingshot", quietly=TRUE))' \
    || die "slingshot not loadable in $TJ_ENV"
fi

# ---- sgscr_scenic: Python + pySCENIC ----
if ! env_exists "$SC_ENV"; then
  info "creating $SC_ENV (python + pySCENIC)..."
  "$INSTALLER" create -y -n "$SC_ENV" -c conda-forge "python=3.10" numpy pandas cython
  "$INSTALLER" run -n "$SC_ENV" pip install pyscenic loompy
else
  info "$SC_ENV exists; verifying pyscenic importable"
  "$INSTALLER" run -n "$SC_ENV" python -c "import pyscenic" || die "pyscenic not importable in $SC_ENV"
fi

info "all three isolated envs ready: $CC_ENV / $TJ_ENV / $SC_ENV"
