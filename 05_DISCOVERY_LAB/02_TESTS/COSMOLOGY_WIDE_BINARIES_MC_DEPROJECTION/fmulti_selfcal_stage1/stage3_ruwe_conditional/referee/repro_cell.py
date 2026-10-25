# Referee: reproducao independente de uma celula (config, rep) do synthetic_identifiability.py
# Uso: python -I -B repro_cell.py CFG REP OUT.json   (apenas sintetico; nenhum dado real)
import sys, json
sys.dont_write_bytecode = True
S3 = "/home/user/TamesisTheoryCompleteResearchArchive/05_DISCOVERY_LAB/02_TESTS/COSMOLOGY_WIDE_BINARIES_MC_DEPROJECTION/fmulti_selfcal_stage1/stage3_ruwe_conditional"
sys.path.insert(0, S3)
import synthetic_identifiability as si
cfg_id, rep, out = int(sys.argv[1]), int(sys.argv[2]), sys.argv[3]
r = si.run_replicate(dict(config_id=cfg_id, rep=rep, cfg=si.CONFIGS[cfg_id], n_mc=si.N_MC, n_boot=si.N_BOOT))
json.dump(r, open(out, "w"))
print("done", cfg_id, rep, r["scalar"]["a0"], r["cond_joint"]["a0"])
