# Referee: contrafactual (IN-SAMPLE, circular, so' ilustrativo): a mesma regra com correcao de vies UNILATERAL (IC/F) em vez de bilateral+L4.
import json, math, sys
import numpy as np
sys.dont_write_bytecode = True
S3 = "/home/user/TamesisTheoryCompleteResearchArchive/05_DISCOVERY_LAB/02_TESTS/COSMOLOGY_WIDE_BINARIES_MC_DEPROJECTION/fmulti_selfcal_stage1/stage3_ruwe_conditional/"
d = json.load(open(S3 + "synthetic_identifiability.json")); R = d["replicates"]; C = d["meta"]["configs"]
A, B = 1.082288e-10, 6.800218e-10; F = d["operating_characteristics"]["F3_bilateral_bias_factor"]; thr = d["cross_config_gates"]["thresholds_95th_pct_of_null"]
def v(l, h):
    a, b = l <= A <= h, l <= B <= h
    return "BOTH_SURVIVE" if a and b else "A_ONLY" if a else "B_ONLY" if b else "BOTH_FALSIFIED"
tot = {}
for c in sorted(R, key=int):
    cnt = {}
    for r in R[c]:
        rc = r["cond_joint"]["recal"]; zmax = float(np.max(np.abs(rc["z_bins0_3"])))
        if rc["T_univ"] > thr["T_univ"] or zmax > thr["zmax"] or rc.get("gof_chi2", 0) > thr["gof_chi2"]: k = "ND[GATE]"
        else:
            cs = r["cond_strata"]; a0 = r["cond_joint"]["a0"]
            if abs(math.log10(a0 / A)) > 1 and abs(math.log10(a0 / B)) > 1: k = "ND[G1]"
            elif (not cs["1"]["converged"]) or cs["1"]["boundary_zero"] or cs["1"]["boundary_top"] or cs["0"]["boundary_top"]: k = "ND[G2]"
            else: k = v(rc["ci_lo"] / F, rc["ci_hi"] / F)   # correcao unilateral (Estagio 2 dividia por 1,4-1,6)
        cnt[k] = cnt.get(k, 0) + 1; tot[k] = tot.get(k, 0) + 1
    print(f"cfg{c:>2} a0_true={'A' if C[c]['a0']==A else 'B' if C[c]['a0']==B else 'Newt'} {C[c]['name']:<42} {cnt}")
print("TOTAL", tot)
