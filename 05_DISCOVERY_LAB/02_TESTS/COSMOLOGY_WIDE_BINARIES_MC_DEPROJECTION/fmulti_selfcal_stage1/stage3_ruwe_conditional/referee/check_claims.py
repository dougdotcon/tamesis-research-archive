# Referee: verificacoes sobre o JSON sintetico (somente leitura; nenhum dado real).
# Uso: python -I -B check_claims.py [json]
import json, math, sys, itertools
import numpy as np
sys.dont_write_bytecode = True
S3 = "/home/user/TamesisTheoryCompleteResearchArchive/05_DISCOVERY_LAB/02_TESTS/COSMOLOGY_WIDE_BINARIES_MC_DEPROJECTION/fmulti_selfcal_stage1/stage3_ruwe_conditional/"
d = json.load(open(sys.argv[1] if len(sys.argv) > 1 else S3 + "synthetic_identifiability.json"))
R = d["replicates"]; C = d["meta"]["configs"]
A, B = 1.082288e-10, 6.800218e-10
LOGAB = math.log10(B / A)
F = d["operating_characteristics"]["F3_bilateral_bias_factor"]
thr = d["cross_config_gates"]["thresholds_95th_pct_of_null"]
print("5784/30203 =", 5784 / 30203, "| log10(aB/aA) =", LOGAB, "| 2log10F =", 2 * math.log10(F), "| w_max p/ L4 =", LOGAB - 2 * math.log10(F))

def rule(r, thr=thr, F=F, use_L4=True, use_gates=True):
    rc = r["cond_joint"]["recal"]; lo, hi = rc["ci_lo"], rc["ci_hi"]; a0 = r["cond_joint"]["a0"]
    if lo is None or a0 is None: return "NO_FIT"
    zmax = float(np.max(np.abs(rc["z_bins0_3"])))
    if use_gates and (rc["T_univ"] > thr["T_univ"] or zmax > thr["zmax"] or rc.get("gof_chi2", 0) > thr["gof_chi2"]): return "ND[GATE]"
    cs = r["cond_strata"]
    if abs(math.log10(a0 / A)) > 1 and abs(math.log10(a0 / B)) > 1: return "ND[G1]"
    if (not cs["1"]["converged"]) or cs["1"]["boundary_zero"] or cs["1"]["boundary_top"] or cs["0"]["boundary_top"]: return "ND[G2]"
    cut = r["cut_low_only"]["a0"]
    if cut is not None and abs(math.log10(a0 / cut)) > 0.3: return "ND[G9]"
    def v(l, h):
        a, b = l <= A <= h, l <= B <= h
        return "BOTH_SURVIVE" if a and b else "A_ONLY" if a else "B_ONLY" if b else "BOTH_FALSIFIED"
    m = v(lo, hi)
    if m == "BOTH_SURVIVE": return "ND[CI]"
    if use_L4 and math.log10(hi * F / (lo / F)) >= LOGAB: return "ND[LARG]"
    if v(lo / F, hi * F) != m: return "ND[VIES]"
    return m

tot = {}
for c in sorted(R, key=int):
    cnt = {}
    for r in R[c]:
        k = rule(r); cnt[k] = cnt.get(k, 0) + 1; tot[k] = tot.get(k, 0) + 1
    print("cfg", c, C[c]["name"], cnt)
print("TOTAL regra completa (140 reps):", tot)

# 1) impossibilidade estrutural de A_ONLY: lo<=A exige meia-largura>=vies; L4 exige largura<0,363
print("\n== A_ONLY alcancavel? (H_A verdadeiro): lo<=a0_A e largura<%.3f dex" % (LOGAB - 2 * math.log10(F)))
for c in sorted(R, key=int):
    if C[c]["a0"] != A: continue
    n_lo = n_w = n_both = 0; ws = []; bs = []
    for r in R[c]:
        rc = r["cond_joint"]["recal"]; lo, hi = rc["ci_lo"], rc["ci_hi"]
        if lo is None: continue
        w = math.log10(hi / lo); ws.append(w); bs.append(math.log10(r["cond_joint"]["a0"] / A))
        a = lo <= A; b = w < LOGAB - 2 * math.log10(F)
        n_lo += a; n_w += b; n_both += (a and b)
    print(f"cfg{c:>2}: lo<=aA {n_lo}/10 | largura<0.363 {n_w}/10 | ambos {n_both}/10 | vies medio {np.mean(bs):+.3f} | largura media {np.mean(ws):.2f}")

# 1b) condicao analitica: IC aprox. simetrico em log, vies b, meia-largura h: lo<=aA => h>=b ; L4 => 2h < LOGAB-2log10F
b_max = (LOGAB / 4)
print("vies maximo tolerado para A_ONLY ser alcancavel com F=10^b: b < LOGAB/4 = %.4f dex; F^: log10F=%.4f; vies medio observado cfg0,1,12 ~ +0.21" % (b_max, math.log10(F)))

# 2) regra nunca falsifica nada (nem Newtoniano)
print("\nBOTH_FALSIFIED na regra completa em qualquer cfg:", tot.get("BOTH_FALSIFIED", 0), "| A_ONLY:", tot.get("A_ONLY", 0), "| B_ONLY:", tot.get("B_ONLY", 0))
# sem L4/L5 o que acontece
tot2 = {}
for c in R:
    for r in R[c]:
        k = rule(r, use_L4=False); tot2[k] = tot2.get(k, 0) + 1
print("sem L4:", tot2)

# 3) limiares: leave-one-config-out (nulos 0,1,12)
print("\n== Limiares in-sample vs leave-one-config-out (nulos 0,1,12)")
def vals(c): 
    out = []
    for r in R[c]:
        rc = r["cond_joint"]["recal"]; out.append((rc["T_univ"], float(np.max(np.abs(rc["z_bins0_3"]))), rc.get("gof_chi2", 0)))
    return np.array(out)
V = {c: vals(c) for c in ("0", "1", "12")}
allv = np.vstack(list(V.values()))
th = np.percentile(allv, 95, axis=0); print("limiares in-sample (T,zmax,chi2):", th, "vs json", [thr["T_univ"], thr["zmax"], thr["gof_chi2"]])
print("fracao de nulos que dispara QUALQUER gate in-sample:", np.mean(np.any(allv > th, axis=1)))
for out in ("0", "1", "12"):
    tr = np.vstack([V[c] for c in V if c != out]); t = np.percentile(tr, 95, axis=0)
    print(f" holdout cfg{out}: limiares de 20 reps {np.round(t,2)} -> taxa de gate-any em cfg{out}: {np.mean(np.any(V[out] > t, axis=1)):.2f}")
# FP 'any gate' sobre os 3 nulos (in-sample) e Wilson
k = int(np.sum(np.any(allv > th, axis=1))); n = len(allv); p = k / n; z = 1.96
wl = (p + z*z/(2*n) - z*math.sqrt(p*(1-p)/n + z*z/(4*n*n))) / (1 + z*z/n); wu = (p + z*z/(2*n) + z*math.sqrt(p*(1-p)/n + z*z/(4*n*n))) / (1 + z*z/n)
print(f"FP any-gate in-sample {k}/{n}; Wilson [{wl:.2f},{wu:.2f}]")
# bootstrap: incerteza do percentil 95 com n=30
rng = np.random.default_rng(1); bt = np.array([np.percentile(allv[rng.integers(0, 30, 30)], 95, axis=0) for _ in range(2000)])
print("IC90 bootstrap do limiar (T,zmax,chi2): lo", np.round(np.percentile(bt, 5, axis=0), 2), "hi", np.round(np.percentile(bt, 95, axis=0), 2))

# 4) vies F3 e sua origem
bs = []
for c in ("0", "1", "12"):
    bs.append(np.median([math.log10(r["cond_joint"]["a0"] / A) for r in R[c]]))
print("\nvies mediano cond por cfg 0,1,12:", np.round(bs, 4), "-> max=", max(abs(x) for x in bs), "10^max =", 10 ** max(abs(x) for x in bs))
allb = [math.log10(r["cond_joint"]["a0"] / A) for c in ("0", "1", "12") for r in R[c]]
print("vies individual: max", max(allb), "frac reps com vies> log10F:", np.mean(np.array(allb) > math.log10(F)))

# 5) contagens Estagio 2 (escalar + camada)
print("\n== Estagio 2 (stage2_scalar) por cfg (JSON oc):")
for c, v in d["operating_characteristics"]["by_config"].items():
    print(c, v["name"], v["counts"]["stage2_scalar"])
# 6) cfg9 bias
print("\ncfg9 vies cond:", np.round([math.log10(r["cond_joint"]["a0"] / B) for r in R["9"]], 3))
print("cfg9 CI recal (lo,hi)/aB:", [(round(r["cond_joint"]["recal"]["ci_lo"] / B, 2), round(r["cond_joint"]["recal"]["ci_hi"] / B, 2)) for r in R["9"]])
