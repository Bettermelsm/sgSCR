#!/usr/bin/env python
# s9_scenic.py — SCENIC 调控子分析（需求 15，隔离环境 sgscr_scenic，pySCENIC 三步）
# Validates: Requirements 15.1-15.2 (P11)
#
# 前置校验（P11）: 物种 ranking feather 存在、TF 清单存在、主环境 R 可用于导出。
# 流程: 主环境导出表达矩阵 → pyscenic grn → ctx → aucell →
#       scenic/regulon_activity.csv + regulon_<TF>.png（top 调控子热图）
# 用法: python s9_scenic.py --config <yaml>

import argparse
import csv
import os
import shutil
import subprocess
import sys
from datetime import datetime

import pandas as pd
import yaml

STEP = "s9_scenic"
STARTED_AT = None   # main() 开头赋值，manifest 条目用


def manifest_append(cfg, status, detail, outputs, params):
    """manifest 条目追加（Python 侧与 R 侧 sgs_manifest_append_step 字段同构，只追加）"""
    import json
    from datetime import datetime
    p = cfg["paths"]["manifest_path"]
    m = json.load(open(p, encoding="utf-8"))
    finished = datetime.now().astimezone().isoformat(timespec="seconds")
    dur = (datetime.fromisoformat(finished) - datetime.fromisoformat(STARTED_AT)).total_seconds()
    try:
        from importlib.metadata import version as _v
        pv = {"pyscenic": _v("pyscenic")}
    except Exception:
        pv = {}
    m["steps"].append({
        "step_id": STEP, "started_at": STARTED_AT, "finished_at": finished,
        "duration_sec": dur, "env_name": cfg["envs"]["scenic"],
        "package_versions": pv,
        "quality_gate": {"status": status, "detail": detail},
        "human_confirmation": None, "resolved_params": params, "outputs": outputs,
    })
    json.dump(m, open(p, "w", encoding="utf-8"), ensure_ascii=False, indent=2)


def prefail(missing):
    if missing:
        print(f"[{STEP}] 前置校验失败，缺失必需项: {', '.join(sorted(missing))}")
        sys.exit(1)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--config", required=True)
    args = ap.parse_args()
    cfg = yaml.safe_load(open(args.config, encoding="utf-8"))
    global STARTED_AT
    STARTED_AT = datetime.now().astimezone().isoformat(timespec="seconds")

    feather = os.path.expanduser(cfg["scenic"]["species_db_path"])
    tf_file = os.path.expanduser(cfg["scenic"].get("tf_list_path", os.path.join(
        os.path.dirname(feather), "allTFs_hg38.txt")))
    motif_tbl = os.path.expanduser(cfg["scenic"].get("motif_annotations_path", os.path.join(
        os.path.dirname(feather), "motifs-v9-nr.hgnc-m0.001-o0.0.tbl")))
    workdir = os.path.join(os.path.dirname(cfg["paths"]["manifest_path"]), "..", "scenic")
    problems = []
    if not os.path.isfile(feather):
        problems.append(f"物种参考数据库 feather: {feather}")
    if not os.path.isfile(tf_file):
        problems.append(f"TF 清单: {tf_file}")
    if not os.path.isfile(motif_tbl):
        problems.append(f"motif 注释库: {motif_tbl}")
    if shutil.which("Rscript") is None and not os.path.exists(
            os.path.expanduser(f"~/anaconda3/envs/{cfg['envs']['main']}/bin/Rscript")):
        problems.append(f"主环境Rscript({cfg['envs']['main']})")
    prefail(problems)

    os.makedirs(workdir, exist_ok=True)
    conda_rscript = os.path.expanduser(
        f"~/anaconda3/envs/{cfg['envs']['main']}/bin/Rscript")
    main_rscript = conda_rscript if os.path.exists(conda_rscript) else shutil.which("Rscript")
    skill_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

    # ---- 1. 主环境导出表达矩阵 ----
    export = os.path.join(skill_dir, "scripts", "common", "export_scenic_input.R")
    subprocess.run([main_rscript, export, "--config", args.config,
                    "--outdir", workdir], check=True)
    expr_csv = os.path.join(workdir, "expr.csv")
    # ctx 的 modules_from_adjacencies 要求 cells×genes（与 grn/aucell 的 -t 方向相反），
    # 单独落盘转置矩阵；否则要么断言报错、要么静默产出 0 调控子
    expr_t_csv = os.path.join(workdir, "expr_cells.csv")
    if not os.path.exists(expr_t_csv):
        pd.read_csv(expr_csv, index_col=0).T.to_csv(expr_t_csv)

    n_workers = "2"
    adj_csv = os.path.join(workdir, "adjacencies.tsv")
    reg_csv = os.path.join(workdir, "regulons.csv")
    auc_csv = os.path.join(workdir, "auc_matrix.csv")

    # ---- 2-4. pySCENIC 三步（grn / ctx / aucell） ----
    pyscenic = shutil.which("pyscenic") or os.path.join(
        os.path.dirname(sys.executable), "pyscenic")   # env bin 下的 console script
    if not os.path.exists(pyscenic):
        print(f"[{STEP}] 前置校验失败，缺失必需项: pyscenic 可执行文件")
        sys.exit(1)
    def run_scenic(sub, *a):
        cmd = (pyscenic if isinstance(pyscenic, list) else [pyscenic]) + [sub, *a]
        print(f"[{STEP}] {' '.join(cmd)}")
        subprocess.run(cmd, check=True)

    def valid(out):   # 判活按行数：只读到第 2 行即可区分"仅表头"与"有实际数据"
        if not os.path.exists(out) or os.path.getsize(out) == 0:
            return False
        try:
            with open(out, "r", encoding="utf-8", errors="replace") as fh:
                for i, _ in enumerate(fh):
                    if i >= 1:
                        return True
        except OSError:
            return False
        return False

    # 各步产物存在即复用（幂等）：GRNBoost 并行有非确定性，重跑可能得到不同结果
    if valid(adj_csv):
        print(f"[{STEP}] 复用已有产物: {adj_csv}")
    else:
        # 固定 seed 保证 GRNBoost 回归器随机状态可复现
        run_scenic("grn", "-t", expr_csv, tf_file, "-o", adj_csv, "--num_workers",
                   n_workers, "--seed", str(cfg["run"]["random_seed"]), "--client", "local")
    # ctx: .tsv adjacency 走 adjacencies2modules 分支，必须显式提供 cells×genes 表达矩阵
    if valid(reg_csv):
        print(f"[{STEP}] 复用已有产物: {reg_csv}")
    else:
        run_scenic("ctx", adj_csv, feather, "--annotations_fname", motif_tbl,
                   "--expression_mtx_fname", expr_t_csv, "-o", reg_csv,
                   "--nes_threshold", str(cfg["scenic"].get("nes_threshold", 3.0)),
                   "--num_workers", n_workers, "--client", "local")
    # 清洗：TargetGenes 为空的签名会让 aucell 报 "must have at least one gene"
    reg = pd.read_csv(reg_csv, header=[0, 1], index_col=[0, 1])
    tg = reg[("Enrichment", "TargetGenes")].astype(str)
    bad = reg[tg.isin(["nan", "", "[]", "()"])].index
    if len(bad):
        print(f"[{STEP}] 移除 TargetGenes 为空的签名 {len(bad)} 条: {list(bad)[:3]}...")
        reg.drop(index=bad).to_csv(reg_csv)
    # aucell 子命令不接受 --client（grn/ctx 专属）
    if valid(auc_csv):
        print(f"[{STEP}] 复用已有产物: {auc_csv}")
    else:
        run_scenic("aucell", "-t", expr_csv, reg_csv, "-o", auc_csv,
                   "--num_workers", n_workers)

    # ---- 5. regulon 活性落盘 + top 调控子热图 ----
    # aucell 输出 regulons×cells，统一转置为 cells×regulons（regulon_activity.csv 落盘此方向）
    auc = pd.read_csv(auc_csv, index_col=0).T
    auc.to_csv(os.path.join(workdir, "regulon_activity.csv"))
    if auc.shape[1] == 0:
        print(f"[{STEP}] Quality_Gate: FAIL - 未检出任何调控子")
        manifest_append(cfg, "FAIL", "未检出任何调控子",
                        outputs=[os.path.join(workdir, "regulon_activity.csv")], params={})
        sys.exit(1)

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    top = auc.mean().sort_values(ascending=False).head(min(20, auc.shape[1])).index
    fig, ax = plt.subplots(figsize=(cfg["figure"]["width_in"] + 4,
                                    cfg["figure"]["height_in"] + 2))
    im = ax.imshow(auc[top].T.values, aspect="auto", cmap="viridis")
    ax.set_yticks(range(len(top))); ax.set_yticklabels(top, fontsize=6)
    ax.set_title(f"Regulon activity (top {len(top)}, n={auc.shape[0]} cells)")
    fig.colorbar(im, ax=ax, label="AUC")
    fig.tight_layout()
    fig.savefig(os.path.join(workdir, "regulon_top_heatmap.png"), dpi=cfg["figure"]["dpi"])
    for tf in list(top[:3]):
        fig2, ax2 = plt.subplots(figsize=(cfg["figure"]["width_in"], cfg["figure"]["height_in"]))
        sc = ax2.imshow(auc[tf].values[None, :], aspect="auto", cmap="viridis")
        ax2.set_yticks([0]); ax2.set_yticklabels([tf])
        ax2.set_xticks([])
        fig2.colorbar(sc, ax=ax2, label="AUC")
        safe = tf.replace("/", "_").replace("(", "").replace(")", "").replace("+", "p").replace("-", "m")
        fig2.savefig(os.path.join(workdir, f"regulon_{safe}.png"), dpi=cfg["figure"]["dpi"])
        plt.close(fig2)

    print(f"[{STEP}] Quality_Gate: PASS - {auc.shape[1]} regulons x {auc.shape[0]} cells")
    manifest_append(cfg, "PASS", f"{auc.shape[1]} regulons x {auc.shape[0]} cells",
                    outputs=[os.path.join(workdir, "regulon_activity.csv"),
                             os.path.join(workdir, "regulon_top_heatmap.png")],
                    params={"nes_threshold": cfg["scenic"].get("nes_threshold", 3.0),
                            "target_subset": cfg["scenic"].get("target_subset"),
                            "n_regulons": int(auc.shape[1])})
    print(f"[{STEP}] S9 complete")


if __name__ == "__main__":
    main()
