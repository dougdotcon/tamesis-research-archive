# Referee: sensibilidade do poder/FP dos gates a incerteza do limiar (30 nulos), so' a partir do JSON.
import json, sys, math
import numpy as np
sys.dont_write_bytecode = True
S3 = "/home/user/TamesisTheoryCompleteResearchArchive/05_DISCOVERY_LAB/02_TESTS/COSMOLOGY_WIDE_BINARIES_MC_DEPROJECTION/fmulti_selfcal_stage1/stage3_ruwe_conditional/"
d = json.load(open(S3 + "synthetic_identifiability.json")); R = d["replicates"]
def vals(c): return np.array([[r["cond_joint"]["recal"]["T_univ"], float(np.max(np.abs(r["cond_joint"]["recal"]["z_bins0_3"]))), r["cond_joint"]["recal"].get("gof_chi2", 0)] for r in R[c]])
null = np.vstack([vals(c) for c in ("0", "1", "12")])
rng = np.random.default_rng(7)
B = np.array([np.percentile(null[rng.integers(0, 30, 30)], 95, axis=0) for _ in range(4000)])
th0 = np.percentile(null, 95, axis=0)
def anyg(v, th): return np.any(v > th, axis=1).mean()
print("limiar in-sample", np.round(th0, 2))
for c in ("3", "4", "13", "5", "8", "10", "11"):
    v = vals(c)
    pw = np.array([anyg(v, th) for th in B])
    print(f"cfg{c:>2} {d['meta']['configs'][c]['name']:<42} taxa 'qualquer gate' in-sample {anyg(v, th0):.1f} | sob incerteza do limiar: mediana {np.median(pw):.1f}, IC90 [{np.percentile(pw,5):.1f},{np.percentile(pw,95):.1f}]")
fp = np.array([anyg(null, th) for th in B]); print("FP qualquer-gate nos 30 nulos: in-sample", anyg(null, th0), "| (otimismo: limiar e FP medidos nos mesmos dados)")
# FP esperado fora da amostra: reamostra 30 nulos para limiar e avalia em reamostra independente de 30 (bootstrap cruzado)
rng = np.random.default_rng(11); oos = []
for _ in range(4000):
    tr = null[rng.integers(0, 30, 30)]; te = null[rng.integers(0, 30, 30)]
    oos.append(anyg(te, np.percentile(tr, 95, axis=0)))
print("FP qualquer-gate fora-da-amostra (bootstrap cruzado): media", np.mean(oos), "IC90", np.percentile(oos, [5, 95]))
