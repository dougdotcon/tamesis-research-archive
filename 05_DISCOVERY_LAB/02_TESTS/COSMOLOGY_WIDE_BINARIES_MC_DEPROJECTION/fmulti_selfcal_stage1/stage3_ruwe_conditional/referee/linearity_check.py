# Referee: a aproximacao LINEAR em f do FPB (U3 validou so' t=0,47 interior) e' boa em EXTRAPOLACAO (t ate 3)?
# Estrato RUWE-alto sintetico (f_true=0,45), N=12944 -> ~2480 sistemas; seeds 20260955900-905 (bloco reservado, nao usadas pelo rascunho).
import sys, numpy as np
sys.dont_write_bytecode = True
S3 = "/home/user/TamesisTheoryCompleteResearchArchive/05_DISCOVERY_LAB/02_TESTS/COSMOLOGY_WIDE_BINARIES_MC_DEPROJECTION/fmulti_selfcal_stage1/stage3_ruwe_conditional"
sys.path.insert(0, S3)
import synthetic_identifiability as si
rng = np.random.default_rng(20260955900)
cat = si.generate_equal_bin_catalog(12944, rng)
n = cat["n_sys"]; rt = np.random.default_rng(20260955901)
label = (rt.random(n) < si.PI_HIGH).astype(int)
f_true = np.where(label == 1, 0.45, 0.02)
tr = si.build_truth(cat, f_true, si.A0_A, rt)
idx = np.where(label == 1)[0]; ck = si.subset_cat(cat, idx); v = tr["v_p"][idx]
def delta(f): return np.array(si.run_pipeline(ck, v, f, 20260955902, 60)["delta_obs_newt_primary"])
fa, fb = 0.40, 0.55   # par pareado (passo 0,15 como no FPB)
da, db = delta(fa), delta(fb)
print("n estrato alto:", len(idx), "| delta(fa) =", np.round(da, 3), "| delta(fb) =", np.round(db, 3))
for f in (0.45, 0.30, 0.20, 0.10, 0.70, 0.90):
    t = (f - fa) / (fb - fa); lin = da + t * (db - da); ex = delta(f)
    print(f"f={f:.2f} t={t:+.2f} erro linear por bin (lin-exato): {np.round(lin - ex, 3)} | max|erro|={np.max(np.abs(lin-ex)):.3f} dex")
