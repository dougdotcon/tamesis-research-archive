# Referee: compara replicas reproduzidas com o JSON registrado (diferenca absoluta/relativa maxima, campo a campo).
import json, sys
sys.dont_write_bytecode = True
S3 = "/home/user/TamesisTheoryCompleteResearchArchive/05_DISCOVERY_LAB/02_TESTS/COSMOLOGY_WIDE_BINARIES_MC_DEPROJECTION/fmulti_selfcal_stage1/stage3_ruwe_conditional/"
ref = json.load(open(S3 + "synthetic_identifiability.json"))["replicates"]
def walk(a, b, path, out):
    if isinstance(a, dict):
        for k in a:
            if k == "elapsed_s": continue
            walk(a[k], b[k], path + "/" + str(k), out)
    elif isinstance(a, list):
        for i, (x, y) in enumerate(zip(a, b)): walk(x, y, path + f"[{i}]", out)
        if len(a) != len(b): out.append((path + " LEN", len(a), len(b), 1e9))
    elif isinstance(a, (int, float)) and not isinstance(a, bool) and b is not None:
        out.append((path, a, b, abs(a - b) / (abs(a) + 1e-300)))
    elif a != b: out.append((path, a, b, 0.0 if a == b else 1e9))
for f in sys.argv[1:]:
    cfg, rep = f.split("/")[-1][1:-5].split("_r")
    new = json.load(open(f)); old = ref[cfg][int(rep)]
    out = []; walk(old, new, "", out)
    mx = max(out, key=lambda t: t[3])
    nz = [t for t in out if t[3] > 0]
    print(f"cfg{cfg} rep{rep}: {len(out)} campos comparados; {len(nz)} com diferenca != 0; max rel diff = {mx[3]:.3e} em {mx[0]}")
    for t in sorted(nz, key=lambda t: -t[3])[:5]: print("   ", t)
    for key in ("scalar/a0", "cond_joint/a0", "cond_joint/recal/ci_lo", "cond_joint/recal/ci_hi", "cond_joint/recal/T_univ", "cond_joint/recal/gof_chi2", "cond_strata/0/f_hat", "cond_strata/1/f_hat"):
        o, n = old, new
        for p in key.split("/"): o, n = o[p], n[p]
        print(f"   {key:<28} registrado {o!r:<24} reproduzido {n!r}")
