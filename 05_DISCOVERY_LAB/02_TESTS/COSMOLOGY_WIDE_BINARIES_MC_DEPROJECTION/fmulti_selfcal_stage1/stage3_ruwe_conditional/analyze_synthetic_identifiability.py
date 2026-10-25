#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Tabelas derivadas de synthetic_identifiability.json (somente leitura do JSON sintetico).

Nenhum dado real. Nenhuma simulacao nova: apenas agrega o que synthetic_identifiability.py gravou
e calcula comparacoes entre configs (testes de duas amostras sobre o vies em dex).
Todo numero citado no PREREGISTRATION_STAGE3_DRAFT.md que nao esteja em synthetic_identifiability.log
vem de analysis_summary.log (saida deste script).

Uso: python -I -B analyze_synthetic_identifiability.py [--json caminho] [--out caminho]
"""
from __future__ import annotations

import argparse
import json
import math
import sys
from pathlib import Path

sys.dont_write_bytecode = True
import numpy as np  # noqa: E402
from scipy import stats  # noqa: E402

HERE = Path(__file__).resolve().parent
A0_A = 1.082288e-10
A0_B = 6.800218e-10


def lb(a0, a0_true):
    if a0 is None or a0_true is None or a0 <= 0:
        return None
    return math.log10(a0 / a0_true)


def arr(xs):
    return np.array([x for x in xs if x is not None and np.isfinite(x)], float)


def fmt(a, nd=3):
    if len(a) == 0:
        return "n/a"
    sd = a.std(ddof=1) if len(a) > 1 else float("nan")
    return f"{a.mean():+.{nd}f} +-{sd:.{nd}f} (med {np.median(a):+.{nd}f}, n={len(a)})"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--json", default=str(HERE / "synthetic_identifiability.json"))
    ap.add_argument("--out", default=str(HERE / "analysis_summary.log"))
    args = ap.parse_args()
    d = json.load(open(args.json, encoding="utf-8"))
    raw = d["replicates"]
    cfgs = d["meta"]["configs"]
    outf = open(args.out, "w", encoding="utf-8")

    def log(m=""):
        print(m, flush=True)
        outf.write(m + "\n")
        outf.flush()

    log("analyze_synthetic_identifiability.py -- tabelas derivadas (somente leitura de synthetic_identifiability.json)")
    log(f"json: {args.json}")
    log(f"configs presentes: {sorted(raw, key=int)}; reps por config: {[len(raw[c]) for c in sorted(raw, key=int)]}")
    log("")

    # ------------------------------------------------------------------ tabela 1: vies por metodo
    log("=== T1. Vies de log10(a0_fit/a0_true) por metodo (dex): media +-sd (mediana, n) ===")
    for c in sorted(raw, key=int):
        cf = cfgs[c]
        a0t = cf["a0"]
        rs = raw[c]
        if a0t is None:
            continue
        s = arr([lb(r["scalar"]["a0"], a0t) for r in rs])
        j = arr([lb(r["cond_joint"]["a0"], a0t) for r in rs])
        l = arr([lb(r["cut_low_only"]["a0"], a0t) for r in rs])
        log(f"cfg{c:>2} {cf['name']:<40} N={cf['n_sys']:>5} | escalar {fmt(s)} | cond-conjunto {fmt(j)} | corte-baixo {fmt(l)}")
    log("")

    # ------------------------------------------------------------------ tabela 2: largura de IC e cobertura
    log("=== T2. Largura do IC95% (log10 hi/lo) e cobertura do a0 verdadeiro ===")
    for c in sorted(raw, key=int):
        cf = cfgs[c]
        a0t = cf["a0"]
        rs = raw[c]

        def width(lo, hi):
            return math.log10(hi / lo) if (lo and hi and lo > 0) else None

        def cov(lo, hi):
            return (lo <= a0t <= hi) if (a0t and lo and hi) else None
        ws = arr([width(r["scalar"]["ci_lo"], r["scalar"]["ci_hi"]) for r in rs])
        wr = arr([width(r["cond_joint"]["recal"]["ci_lo"], r["cond_joint"]["recal"]["ci_hi"]) for r in rs])
        wf = arr([width(r["cond_joint"]["fixed"]["ci_lo"], r["cond_joint"]["fixed"]["ci_hi"]) for r in rs])
        wl = arr([width(r["cut_low_only"]["ci_lo"], r["cut_low_only"]["ci_hi"]) for r in rs])
        covs = [cov(r["scalar"]["ci_lo"], r["scalar"]["ci_hi"]) for r in rs]
        covr = [cov(r["cond_joint"]["recal"]["ci_lo"], r["cond_joint"]["recal"]["ci_hi"]) for r in rs]
        covf = [cov(r["cond_joint"]["fixed"]["ci_lo"], r["cond_joint"]["fixed"]["ci_hi"]) for r in rs]
        covl = [cov(r["cut_low_only"]["ci_lo"], r["cut_low_only"]["ci_hi"]) for r in rs]

        def cr(x):
            x = [v for v in x if v is not None]
            return f"{sum(x)}/{len(x)}" if x else "n/a"
        log(f"cfg{c:>2} {cf['name']:<40} | largura(dex) escalar {ws.mean():.2f} cond-recal {wr.mean():.2f} cond-fixo {wf.mean():.2f} corte-baixo {wl.mean():.2f}"
            f" | cobertura escalar {cr(covs)} cond-recal {cr(covr)} cond-fixo {cr(covf)} corte-baixo {cr(covl)}")
    log("")

    # ------------------------------------------------------------------ tabela 3: heterogeneidade e deteccao
    log("=== T3. Deteccao de heterogeneidade (D do escalar; Delta-f do condicional) e f_hat por estrato ===")
    for c in sorted(raw, key=int):
        cf = cfgs[c]
        rs = raw[c]
        n = len(rs)
        kD = sum(r["D_test"]["detects_hetero_95"] for r in rs)
        kF = sum(bool(r["cond_joint"]["recal"]["df_detects_hetero_95"]) for r in rs)
        D = arr([r["D_test"]["D"] for r in rs])
        fl = arr([r["cond_strata"]["0"]["f_hat"] for r in rs])
        fh = arr([r["cond_strata"]["1"]["f_hat"] for r in rs])
        na = arr([r["cond_strata"]["1"]["n_anchor"] for r in rs])
        nl = arr([r["cond_strata"]["0"]["n_anchor"] for r in rs])
        b0 = np.mean([r["cond_strata"]["0"]["boundary_zero"] for r in rs])
        log(f"cfg{c:>2} {cf['name']:<40} N={cf['n_sys']:>5} f_true=({cf['f_low']},{cf['f_high']}) | n_ancora(baixo,alto)=({nl.mean():.0f},{na.mean():.0f})"
            f" | f_hat baixo {fl.mean():.3f}+-{fl.std(ddof=1):.3f} alto {fh.mean():.3f}+-{fh.std(ddof=1):.3f} | frac. f_baixo na borda 0: {b0:.2f}"
            f" | D {D.mean():+.3f} deteccao D {kD}/{n}, deteccao Delta-f {kF}/{n}")
    log("")

    # ------------------------------------------------------------------ testes de duas amostras
    log("=== T4. Comparacoes entre configs (vies do a0 do metodo ESCALAR e do CONDICIONAL; dex) ===")

    def two(c1, c2, getter, label):
        if c1 not in raw or c2 not in raw:
            return
        a = arr([getter(cfgs[c1], r) for r in raw[c1]])
        b = arr([getter(cfgs[c2], r) for r in raw[c2]])
        if len(a) < 2 or len(b) < 2:
            return
        t = stats.ttest_ind(a, b, equal_var=False)
        mw = stats.mannwhitneyu(a, b, alternative="two-sided")
        log(f"{label}: cfg{c1} {a.mean():+.3f} vs cfg{c2} {b.mean():+.3f} (dif {a.mean() - b.mean():+.3f} dex) Welch p={t.pvalue:.3g} MWU p={mw.pvalue:.3g}")
    sc = lambda cf, r: lb(r["scalar"]["a0"], cf["a0"])  # noqa: E731
    cj = lambda cf, r: lb(r["cond_joint"]["a0"], cf["a0"])  # noqa: E731
    cl = lambda cf, r: lb(r["cut_low_only"]["a0"], cf["a0"])  # noqa: E731
    two("0", "1", sc, "escalar, homogeneo(0) vs heterogeneo(1)")
    two("0", "1", cj, "cond-conjunto, homogeneo(0) vs heterogeneo(1)")
    two("12", "1", sc, "escalar, SEM contaminacao(12) vs heterogeneo(1)")
    two("12", "0", sc, "escalar, SEM contaminacao(12) vs homogeneo f=0,10(0)")
    two("1", "10", sc, "escalar, het.(1) vs composicao cresce com gN(10)")
    two("1", "11", sc, "escalar, het.(1) vs composicao decresce com gN(11)")
    two("1", "10", cj, "cond-conjunto, het.(1) vs composicao cresce(10)")
    two("1", "11", cj, "cond-conjunto, het.(1) vs composicao decresce(11)")
    two("1", "3", sc, "escalar, het.(1) vs artefato 0,3 km/s(3)")
    two("1", "3", cj, "cond-conjunto, het.(1) vs artefato 0,3 km/s(3)")
    two("1", "3", cl, "corte-baixo, het.(1) vs artefato 0,3 km/s(3)")
    log("")

    # ------------------------------------------------------------------ gates
    log("=== T5. Gates (universalidade, GOF) e Delta_b medio (alto - baixo, bins 0-3) ===")
    for c in sorted(raw, key=int):
        cf = cfgs[c]
        rs = raw[c]
        Tm = arr([r["cond_joint"]["recal"]["T_univ"] for r in rs])
        G = arr([r["cond_joint"]["recal"].get("gof_chi2") for r in rs])
        Dl = np.nanmean([r["cond_joint"]["recal"]["Delta_point"] for r in rs], axis=0)
        sdD = np.nanmean([r["cond_joint"]["recal"]["Delta_sd"] for r in rs], axis=0)
        emp = np.std([r["cond_joint"]["recal"]["Delta_point"] for r in rs], axis=0, ddof=1)
        log(f"cfg{c:>2} {cf['name']:<40} | Delta_b0-3 medio {np.round(Dl, 3).tolist()} sd_FPB medio {np.round(sdD, 3).tolist()} sd entre-replicatas {np.round(emp, 3).tolist()} razao FPB/entre {np.round(sdD / emp, 2).tolist()} | T {Tm.mean():.2f} (max {Tm.max():.1f}) GOF {G.mean():.2f} (max {G.max():.1f})")
    log("gates com limiar empirico (nulo = configs 0,1,12):")
    log(json.dumps(d.get("cross_config_gates", {}), indent=1, ensure_ascii=False))
    log("")

    # ------------------------------------------------------------------ operating characteristics
    log("=== T6. Caracteristicas operacionais (contagens de veredito) ===")
    oc = d.get("operating_characteristics", {})
    log(f"F3 (fator de vies bilateral) = {oc.get('F3_bilateral_bias_factor')}; faixa Estagio 2 = {oc.get('stage2_F_range')}")
    for c, v in oc.get("by_config", {}).items():
        log(f"cfg{int(c):>2} {v['name']:<40} a0_true={v['a0_true']} n={v['n']}")
        for rule, cnt in v["counts"].items():
            log(f"      {rule:<20} {cnt}")
    log("")

    # ------------------------------------------------------------------ regra completa (parcial)
    log("=== T8. Regra de decisao do rascunho IMPLEMENTADA MAIS COMPLETAMENTE do que em T6 (L1 gates, L2 gatilhos 1/2/9, L3 CI, L4 largura, L5 vies) ===")
    log("    Nao implementados aqui: L0 (pre-condicoes), gatilho 6 (sensibilidade a B/seed), gatilho 10 (variantes), gatilhos 3/8 (so' exigem debunker);")
    log("    gatilho 7 (heterogeneidade ausente) excluido de proposito: nas configs homogeneas (0, 8, 12) ele dispararia por construcao.")
    F3 = oc.get("F3_bilateral_bias_factor", 1.65)
    thr = d.get("cross_config_gates", {}).get("thresholds_95th_pct_of_null", {})
    LOG_AB = math.log10(A0_B / A0_A)

    def verdict(lo, hi):
        a = lo <= A0_A <= hi
        b = lo <= A0_B <= hi
        return "BOTH_SURVIVE" if (a and b) else ("A_ONLY" if a else ("B_ONLY" if b else "BOTH_FALSIFIED"))

    def full_rule(r):
        rc = r["cond_joint"]["recal"]
        lo, hi = rc["ci_lo"], rc["ci_hi"]
        a0 = r["cond_joint"]["a0"]
        if lo is None or a0 is None:
            return "NO_FIT"
        zmax = float(np.max(np.abs(np.array(rc["z_bins0_3"]))))
        if rc["T_univ"] > thr["T_univ"] or zmax > thr["zmax"] or rc.get("gof_chi2", 0) > thr["gof_chi2"]:
            return "NAO_DISTINGUE[GATE]"
        cs = r["cond_strata"]
        t1 = (abs(math.log10(a0 / A0_A)) > 1) and (abs(math.log10(a0 / A0_B)) > 1)
        hi_b = (not cs["1"]["converged"]) or cs["1"]["boundary_zero"] or cs["1"]["boundary_top"]
        lo_b = cs["0"]["boundary_top"]
        cut = r["cut_low_only"]["a0"]
        t9 = (cut is not None) and abs(math.log10(a0 / cut)) > 0.3
        if t1:
            return "NAO_DISTINGUE[GATILHO1]"
        if hi_b or lo_b:
            return "NAO_DISTINGUE[GATILHO2]"
        if t9:
            return "NAO_DISTINGUE[GATILHO9]"
        v = verdict(lo, hi)
        if v == "BOTH_SURVIVE":
            return "NAO_DISTINGUE[CI]"
        if math.log10((hi * F3) / (lo / F3)) >= LOG_AB:
            return "NAO_DISTINGUE[LARGURA]"
        if verdict(lo / F3, hi * F3) != v:
            return "NAO_DISTINGUE[VIES]"
        return v
    for c in sorted(raw, key=int):
        cnt = {}
        for r in raw[c]:
            k = full_rule(r)
            cnt[k] = cnt.get(k, 0) + 1
        log(f"cfg{c:>2} {cfgs[c]['name']:<40} a0_true={cfgs[c]['a0']} | {cnt}")
    log("")

    # ------------------------------------------------------------------ N-sweep
    log("=== T7. Varredura em N (H1 heterogeneo, a0_true = a0_A): configs 7 (3000), 6 (6000), 5 (12944), 1 (30203) ===")
    for c in ("7", "6", "5", "1"):
        if c not in raw:
            continue
        rs = raw[c]
        cf = cfgs[c]
        na = arr([r["cond_strata"]["1"]["n_anchor"] for r in rs])
        nb0 = arr([r["cond_strata"]["1"]["n_per_bin"][0] for r in rs])
        nl = arr([r["cond_strata"]["0"]["n_anchor"] for r in rs])
        fh = arr([r["cond_strata"]["1"]["f_hat"] for r in rs])
        fl = arr([r["cond_strata"]["0"]["f_hat"] for r in rs])
        kD = sum(r["D_test"]["detects_hetero_95"] for r in rs)
        kF = sum(bool(r["cond_joint"]["recal"]["df_detects_hetero_95"]) for r in rs)
        sj = arr([lb(r["cond_joint"]["a0"], cf["a0"]) for r in rs])
        ss = arr([lb(r["scalar"]["a0"], cf["a0"]) for r in rs])
        wr = arr([math.log10(r["cond_joint"]["recal"]["ci_hi"] / r["cond_joint"]["recal"]["ci_lo"]) for r in rs
                  if r["cond_joint"]["recal"]["ci_lo"]])
        fsd = np.nanmean([r["cond_joint"]["recal"]["f_rep_sd"][1] for r in rs])
        log(f"N={cf['n_sys']:>5} (cfg{c}): n_ancora alto={na.mean():.0f} (min {na.min():.0f}) n_bin0 alto={nb0.mean():.0f} n_ancora baixo={nl.mean():.0f}"
            f" | sd(f_hat alto) entre reps={fh.std(ddof=1):.3f}, sd bootstrap medio={fsd:.3f} | sd(f_hat baixo)={fl.std(ddof=1):.3f}"
            f" | deteccao D {kD}/{len(rs)} Delta-f {kF}/{len(rs)} | vies a0 escalar {ss.mean():+.2f}+-{ss.std(ddof=1):.2f}, cond {sj.mean():+.2f}+-{sj.std(ddof=1):.2f}"
            f" | largura IC cond-recal {wr.mean():.2f} dex (a0_B/a0_A = {math.log10(A0_B / A0_A):.2f} dex)")
    log("")
    log("fim.")
    outf.close()


if __name__ == "__main__":
    main()
