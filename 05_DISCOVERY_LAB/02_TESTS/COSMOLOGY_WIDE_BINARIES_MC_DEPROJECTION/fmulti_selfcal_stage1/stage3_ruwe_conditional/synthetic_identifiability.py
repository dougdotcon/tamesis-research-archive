#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Identificabilidade SINTETICA de f_multi condicional a RUWE -- SPARC-FMULTI-STAGE3 (RASCUNHO).

Autoridade: DISC-DEC-155 (so' rascunho de pre-registro; nenhum dado real).
Bloco de seeds reservado: 20260955000-20260955999 (ver registro de seeds abaixo).

PERGUNTAS (todas respondidas SOMENTE com populacao sintetica):
  Q1. O metodo condicional consegue DISTINGUIR uma populacao com f_multi dependente de RUWE
      (heterogenea) de uma populacao homogenea?
  Q2. O a0 recuperado pelo metodo ESCALAR (Estagio 2, travado) e' enviesado na populacao
      heterogenea, enquanto o CONDICIONAL recupera o a0 injetado?
  Q3. Qual o N minimo (por estrato e por bin-ancora) para isso? Onde NAO distingue?
  Q4. O que acontece quando a contaminacao do estrato RUWE-alto NAO obedece ao modelo
      (mis-especificacao: wobble maior que o modelo; ruido astrometrico subestimado)?

REUSO (sem editar nenhum arquivo existente; apenas import):
  - selfcal_pipeline.py: calibrate_f_multi, run_delta_obs_newt_selfcal, fit_a0, bootstrap_a0_refit
  - build_synthetic_population.py: generate_structural_catalog (catalogo estrutural sintetico)
  - companion_injection.py: sample_companion_injection_mc (aceita f_multi vetorial por broadcast)
  - deprojection_common.py: constantes e amostragem orbital
  - delta_obs_newt.py: apenas a constante BIN_EDGES_LOG_GN_SPARC003 (nenhuma funcao de dado real)
  O gerador de "verdade" (build_truth) replica build_fake_real_dataset com f vetorial
  e dois botoes de mis-especificacao; o teste unitario U1 prova equivalencia exata
  com build_fake_real_dataset no caso escalar (mesma seed => mesmo v_p).

DADO REAL: este script NAO le, abre, nem referencia nenhum arquivo de dado real, nem o
split discovery/holdout. Todos os numeros de entrada vem de (a) constantes citadas com
arquivo:linha ou (b) botoes de cenario sinteticos DECLARADOS como tais abaixo.

NOTA SOBRE SEEDS: as seeds passadas por este script estao todas no bloco reservado.
O pipeline travado (selfcal_pipeline.run_delta_obs_newt_selfcal) deriva, a partir da seed
recebida, sub-geradores por deslocamentos fixos (+111_111_111, +222_222_222, +333_333_333,
+444_444_444, +555_555_555) e calibrate_f_multi usa seed+777 -- esses valores DERIVADOS
caem fora do bloco numerico literal, mas sao funcoes deterministas das seeds registradas.
Isto e' registrado como inconsistencia/ressalva no relatorio (STAGE3_DRAFT_NOTES.md).

Execucao (Python do venv do scratchpad, modo isolado, sem bytecode):
  python -I -B synthetic_identifiability.py --part all --workers 4 --outdir <dir>
"""
from __future__ import annotations

import argparse
import json
import math
import os
import sys
import time
import warnings
from pathlib import Path

os.environ.setdefault("OMP_NUM_THREADS", "1")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
os.environ.setdefault("MKL_NUM_THREADS", "1")
sys.dont_write_bytecode = True  # nao escrever __pycache__ fora do diretorio permitido

HERE = Path(__file__).resolve().parent
STAGE1_DIR = HERE.parent                      # .../fmulti_selfcal_stage1
MC_DIR = STAGE1_DIR.parent                    # .../COSMOLOGY_WIDE_BINARIES_MC_DEPROJECTION
sys.path.insert(0, str(STAGE1_DIR / "analysis"))
sys.path.insert(0, str(MC_DIR / "analysis"))

import numpy as np  # noqa: E402
from scipy import stats  # noqa: E402
from scipy.optimize import curve_fit  # noqa: E402

import deprojection_common as dc  # noqa: E402
import companion_injection as ci  # noqa: E402
import selfcal_pipeline as sp  # noqa: E402
import build_synthetic_population as bp  # noqa: E402
import delta_obs_newt as don  # noqa: E402  (so' a constante de bordas de bin)

warnings.filterwarnings("ignore", category=RuntimeWarning)

# ----------------------------------------------------------------------------------------
# Constantes (cada uma com origem)
# ----------------------------------------------------------------------------------------
BIN_EDGES = np.array(don.BIN_EDGES_LOG_GN_SPARC003, dtype=float)
# PREREGISTRATION_STAGE2.md:124 -- [-11,7012; -9,1728; -8,4667; -7,9752; -7,5548; -6,5224]
assert np.allclose(BIN_EDGES, [-11.7012, -9.1728, -8.4667, -7.9752, -7.5548, -6.5224])
ANCHOR = 4          # PREREGISTRATION_STAGE2.md:134-138 (bin-ancora = 4)
NB = 5
A0_A = 1.082288e-10  # PREREGISTRATION_STAGE2.md:66 (H_A = cH0/2pi)
A0_B = 6.800218e-10  # PREREGISTRATION.md:72; PREREGISTRATION_STAGE2.md:69,319 (H_B = cH0)
N_DISCOVERY = 30203  # PREREGISTRATION_STAGE2.md:87 / RESULTS_PRIMARY_STAGE2.md:25
N_HOLDOUT_DOCUMENTED = 12944  # contagem documentada (PREREGISTRATION.md:100); NENHUM dado do holdout lido
# fracao RUWE-alto (max(RUWE1,RUWE2)>1,4) na amostra de descoberta:
# ADVERSARIAL_DEBUNKER_REPORT.md:179 -> 0,1915 (5784/30203); RESULTS_PRIMARY_STAGE2.md:225-226 -> 19,15%
PI_HIGH = 5784 / 30203
F_HAT_STAGE2 = 0.103666  # RESULTS_PRIMARY_STAGE2.md:78 (so' para comparacao/resemblance)

SEED_BLOCK_LO = 20260955000
SEED_BLOCK_HI = 20260955999
N_REPS_PER_CONFIG = 10
N_CONFIGS = 14
SEEDS_PER_REP = 6  # cat, truth, calib_scalar, calib_low, calib_high, boot
BIN0_MEDIAN_GN_TARGET = 1.8807e-10  # RESULTS_PRIMARY_STAGE2.md:95 (gN mediano do bin 0 REAL, valor ja' publicado)


def seed_for(config_id: int, rep: int, role: int) -> int:
    """seed = 20260955000 + 6*(10*config_id + rep) + role  (config_id<=13 => max 20260955839)."""
    assert 0 <= config_id < N_CONFIGS and 0 <= rep < N_REPS_PER_CONFIG and 0 <= role < SEEDS_PER_REP
    s = SEED_BLOCK_LO + SEEDS_PER_REP * (N_REPS_PER_CONFIG * config_id + rep) + role
    assert SEED_BLOCK_LO <= s <= SEED_BLOCK_HI
    return s


DEV_SEED_BASE = 20260955960  # 960..999 reservadas p/ testes unitarios e calibracao de cenario


def dev_seed(k: int) -> int:
    s = DEV_SEED_BASE + k
    assert 20260955960 <= s <= SEED_BLOCK_HI
    return s


# ----------------------------------------------------------------------------------------
# Cenarios SINTETICOS (botoes declarados; NAO sao afirmacoes sobre dado real)
# ----------------------------------------------------------------------------------------
# f_low/f_high: "f_multi por estrato RUWE" verdadeiro injetado. f_low=0,02 e f_high=0,45 sao
# ext_kms = jitter gaussiano de velocidade (km/s) NAO modelado, adicionado so' a' verdade do estrato RUWE-alto
# (artefato astrometrico hipotetico); 0,3 km/s foi escolhido por exploracao (parte 'explore') por reproduzir
# qualitativamente o perfil do debunker; NAO e' uma medida do dado real.
# valores ESCOLHIDOS para que o analogo sintetico se PAREÇA com o padrao relatado pelo debunker
# (subgrupo RUWE-alto com f na faixa 0,25-0,47 da literatura citada em PREREGISTRATION_STAGE2.md:35,380,
# subgrupo RUWE-baixo quase limpo, prevalencia pooled ~0,10 ~ F_HAT_STAGE2). Isto NAO e' inferencia
# sobre o dado real; a verificacao de semelhanca esta na Parte "calib" (log).
PI_INC = [0.04, 0.08, 0.16, 0.27, 0.41]  # media 0,192 ~ PI_HIGH; fracao RUWE-alto CRESCE com gN (cenario hipotetico H2+)
PI_DEC = PI_INC[::-1]                      # fracao RUWE-alto DECRESCE com gN (cenario hipotetico H2-)
CONFIGS = {
    0: dict(name="H0_homog_f0.10_a0A", f_low=0.10, f_high=0.10, a0=A0_A, n_sys=N_DISCOVERY, wob=1.0, kap=1.0, ext_kms=0.0, pi_bin=None),
    1: dict(name="H1_hetero_a0A", f_low=0.02, f_high=0.45, a0=A0_A, n_sys=N_DISCOVERY, wob=1.0, kap=1.0, ext_kms=0.0, pi_bin=None),
    2: dict(name="H1_hetero_newtonian", f_low=0.02, f_high=0.45, a0=None, n_sys=N_DISCOVERY, wob=1.0, kap=1.0, ext_kms=0.0, pi_bin=None),
    3: dict(name="H1a_hetero_artifact0.3kms_a0A", f_low=0.02, f_high=0.45, a0=A0_A, n_sys=N_DISCOVERY, wob=1.0, kap=1.0, ext_kms=0.3, pi_bin=None),
    4: dict(name="H1a_hetero_artifact0.1kms_a0A", f_low=0.02, f_high=0.45, a0=A0_A, n_sys=N_DISCOVERY, wob=1.0, kap=1.0, ext_kms=0.1, pi_bin=None),
    5: dict(name="H1_hetero_a0A_N12944", f_low=0.02, f_high=0.45, a0=A0_A, n_sys=N_HOLDOUT_DOCUMENTED, wob=1.0, kap=1.0, ext_kms=0.0, pi_bin=None),
    6: dict(name="H1_hetero_a0A_N6000", f_low=0.02, f_high=0.45, a0=A0_A, n_sys=6000, wob=1.0, kap=1.0, ext_kms=0.0, pi_bin=None),
    7: dict(name="H1_hetero_a0A_N3000", f_low=0.02, f_high=0.45, a0=A0_A, n_sys=3000, wob=1.0, kap=1.0, ext_kms=0.0, pi_bin=None),
    8: dict(name="H0_homog_f0.10_a0A_N12944", f_low=0.10, f_high=0.10, a0=A0_A, n_sys=N_HOLDOUT_DOCUMENTED, wob=1.0, kap=1.0, ext_kms=0.0, pi_bin=None),
    9: dict(name="H1_hetero_a0B", f_low=0.02, f_high=0.45, a0=A0_B, n_sys=N_DISCOVERY, wob=1.0, kap=1.0, ext_kms=0.0, pi_bin=None),
    10: dict(name="H2inc_composition_rises_with_gN_a0A", f_low=0.02, f_high=0.45, a0=A0_A, n_sys=N_DISCOVERY, wob=1.0, kap=1.0, ext_kms=0.0, pi_bin=PI_INC),
    11: dict(name="H2dec_composition_falls_with_gN_a0A", f_low=0.02, f_high=0.45, a0=A0_A, n_sys=N_DISCOVERY, wob=1.0, kap=1.0, ext_kms=0.0, pi_bin=PI_DEC),
    12: dict(name="Hclean_f0_a0A", f_low=0.0, f_high=0.0, a0=A0_A, n_sys=N_DISCOVERY, wob=1.0, kap=1.0, ext_kms=0.0, pi_bin=None),
    13: dict(name="H1a_hetero_artifact0.3kms_a0A_N12944", f_low=0.02, f_high=0.45, a0=A0_A, n_sys=N_HOLDOUT_DOCUMENTED, wob=1.0, kap=1.0, ext_kms=0.3, pi_bin=None),
}

N_MC = 100          # (Estagio 1 usou 120; Estagio 2 real usou 200 -- aqui 100 por custo; ver NOTES)
N_BOOT = 300        # bootstrap de IC (Estagio 1: 300-400; Estagio 2 real: 2000)
XTOL = 2e-3         # xtol da bisseccao de f (Estagio 2 real: 5e-4; Estagio 1 validacoes: 5e-3)
F_LO, F_HI = 0.0, 0.9
H_F = 0.15          # passo em f para a interpolacao linear do bootstrap com recalibracao


# ----------------------------------------------------------------------------------------
# Utilitarios de dados sinteticos
# ----------------------------------------------------------------------------------------

def subset_cat(cat: dict, idx: np.ndarray) -> dict:
    out = {}
    for k, v in cat.items():
        if isinstance(v, np.ndarray) and v.shape[:1] == (cat["n_sys"],):
            out[k] = v[idx]
        else:
            out[k] = v
    out["n_sys"] = int(len(idx))
    return out


def generate_equal_bin_catalog(n_sys: int, rng: np.random.Generator, oversample: int = 8) -> dict:
    """Adaptacao declarada: o gerador original (bp.generate_structural_catalog) tem separacao
    log-uniforme, o que NAO reproduz as contagens iguais por bin da amostra real (as bordas de bin
    reais sao quantis iguais; PREREGISTRATION_STAGE2.md:127 -> ~6.040 sistemas por bin).
    Aqui: sorteia `oversample`*n_sys sistemas com o gerador ORIGINAL e retira n_sys/5 de cada bin
    (sem reposicao) pelo gN projetado catalogado. Nenhuma distribuicao nova e' inventada; so' a
    composicao por bin e' igualada. n_sys deve ser multiplo de 5 (ajustado por baixo)."""
    n_sys = (n_sys // NB) * NB
    big = bp.generate_structural_catalog(oversample * n_sys, rng)
    lg = np.log10(sp.G_SI * big["Mtot_cat_kg"] / big["s_m"] ** 2)
    bidx = np.digitize(lg, BIN_EDGES[1:-1])
    inr = (lg >= BIN_EDGES[0]) & (lg <= BIN_EDGES[-1])
    keep = []
    for b in range(NB):
        cand = np.where(inr & (bidx == b))[0]
        if b == 0:
            # O gerador original, com separacao log-uniforme, tem cauda de gN muito baixo no bin 0
            # (mediana sintetica ~5,7e-11 vs 1,8807e-10 REAL, RESULTS_PRIMARY_STAGE2.md:95), o que
            # superestimaria o sinal MOND no bin 0. Aparamos a cauda inferior ate' que a mediana do
            # bin 0 iguale o valor real JA' PUBLICADO (nenhum dado real lido).
            order = cand[np.argsort(lg[cand])]
            lgs = lg[order]
            tgt = np.log10(BIN0_MEDIAN_GN_TARGET)
            L = len(lgs)
            j = 0
            while j < L - n_sys // NB and lgs[j + (L - j) // 2] < tgt:
                j += 1
            cand = order[j:]
        assert len(cand) >= n_sys // NB, (b, len(cand))
        keep.append(rng.choice(cand, size=n_sys // NB, replace=False))
    keep = np.sort(np.concatenate(keep))
    return subset_cat(big, keep)


def build_truth(cat: dict, f_vec, a0_true, rng: np.random.Generator,
                wobble_scale=None, noise_kappa=None, extra_v_si=None) -> dict:
    """Replica bp.build_fake_real_dataset com f_multi VETORIAL (por sistema) e botoes de
    mis-especificacao (escala do wobble e do ruido astrometrico da verdade, por sistema).
    Mesma ordem de chamadas ao rng => equivalencia exata no caso escalar (teste U1)."""
    n_sys = cat["n_sys"]
    s_m = cat["s_m"]
    inj = ci.sample_companion_injection_mc(cat["M1_cat_kg"], cat["M2_cat_kg"], f_vec, n_mc=1, rng=rng)
    mask1, mask2 = inj["mask1"][0], inj["mask2"][0]
    kappa1, kappa2 = inj["kappa_1"][0], inj["kappa_2"][0]
    M1_true, M2_true = inj["M1_true"][0], inj["M2_true"][0]
    Mtot_true = inj["Mtot_true"][0]
    has_multi = inj["has_multi"][0]

    e_sample, _ = dc.sample_eccentricity(cat["e_m"], cat["e_lo"], cat["e_hi"], cat["alpha_ecc"], cat["dpm_sig"], 1, rng)
    i, phi0, phi = dc.sample_orbital_geometry(e_sample, rng)
    v_p_outer_2d, _, _, _ = sp._generate_outer_vp_varying_mass(
        s_m, Mtot_true[None, :], e_sample, i, phi0, phi, a0_boost=a0_true)
    pa = rng.uniform(0.0, 2.0 * np.pi, size=(1, n_sys))
    v_ra = v_p_outer_2d * np.cos(pa)
    v_de = v_p_outer_2d * np.sin(pa)

    ra1, de1 = sp._wobble_component_mc(mask1[None, :], kappa1[None, :], M1_true[None, :], cat["d_mean_pc"], rng)
    ra2, de2 = sp._wobble_component_mc(mask2[None, :], kappa2[None, :], M2_true[None, :], cat["d_mean_pc"], rng)
    if wobble_scale is not None:
        w = np.asarray(wobble_scale)[None, :]
        ra1, de1, ra2, de2 = ra1 * w, de1 * w, ra2 * w, de2 * w
    v_ra = v_ra + ra1 + ra2
    v_de = v_de + de1 + de2
    if extra_v_si is not None:
        # artefato NAO modelado: jitter gaussiano isotropico de velocidade (SI) adicionado so' a' verdade,
        # por sistema (ex.: solucao astrometrica degradada do estrato RUWE-alto). O mock nao o conhece.
        ev = np.asarray(extra_v_si)[None, :]
        v_ra = v_ra + rng.normal(0.0, 1.0, size=v_ra.shape) * ev
        v_de = v_de + rng.normal(0.0, 1.0, size=v_de.shape) * ev
    v_p_clean = np.sqrt(v_ra ** 2 + v_de ** 2)

    sig_ra, sig_de = sp.astrometric_noise_sigma_v_si(
        cat["e_pmRA1"], cat["e_pmRA2"], cat["e_pmDE1"], cat["e_pmDE2"], cat["d_mean_pc"])
    if noise_kappa is not None:
        sig_ra = sig_ra * np.asarray(noise_kappa)
        sig_de = sig_de * np.asarray(noise_kappa)
    v_p_noisy_2d = sp._inject_astrometric_noise(v_p_clean, sig_ra[None, :], sig_de[None, :], rng)
    return {"v_p": v_p_noisy_2d[0], "has_multi": has_multi}


class use_f_vector:
    """Context manager: faz run_delta_obs_newt_selfcal (travado, nao editado) sortear a
    injecao de companheira com f_multi VETORIAL (broadcast em ci.sample_companion_injection_mc,
    que ja' faz `u_has < f_multi`). O argumento escalar `f_multi` da chamada e' ignorado."""

    def __init__(self, f_vec):
        self.f_vec = np.asarray(f_vec, dtype=float)
        self._orig = None

    def __enter__(self):
        self._orig = ci.sample_companion_injection_mc
        f_vec = self.f_vec
        orig = self._orig

        def patched(M1, M2, f_multi, n_mc, rng, *a, **k):
            return orig(M1, M2, f_vec, n_mc, rng, *a, **k)
        ci.sample_companion_injection_mc = patched
        return self

    def __exit__(self, *exc):
        ci.sample_companion_injection_mc = self._orig


def run_pipeline(cat_sub: dict, v_p_sub: np.ndarray, f, seed: int, n_mc: int, f_vec=None):
    args = (cat_sub["s_m"], v_p_sub, cat_sub["M1_cat_kg"], cat_sub["M2_cat_kg"],
            cat_sub["e_m"], cat_sub["e_lo"], cat_sub["e_hi"], cat_sub["alpha_ecc"], cat_sub["dpm_sig"],
            cat_sub["d_mean_pc"], cat_sub["e_pmRA1"], cat_sub["e_pmRA2"], cat_sub["e_pmDE1"], cat_sub["e_pmDE2"])
    kw = dict(f_multi=float(f), bin_edges=BIN_EDGES, n_mc=n_mc, seed=seed, include_wobble=True,
              n_bootstrap=None, return_raw=True)
    if f_vec is None:
        return sp.run_delta_obs_newt_selfcal(*args, **kw)
    with use_f_vector(f_vec):
        return sp.run_delta_obs_newt_selfcal(*args, **kw)


def calibrate(cat_sub: dict, v_p_sub: np.ndarray, seed: int, n_mc: int):
    return sp.calibrate_f_multi(
        cat_sub["s_m"], v_p_sub, cat_sub["M1_cat_kg"], cat_sub["M2_cat_kg"],
        cat_sub["e_m"], cat_sub["e_lo"], cat_sub["e_hi"], cat_sub["alpha_ecc"], cat_sub["dpm_sig"],
        cat_sub["d_mean_pc"], cat_sub["e_pmRA1"], cat_sub["e_pmRA2"], cat_sub["e_pmDE1"], cat_sub["e_pmDE2"],
        bin_edges=BIN_EDGES, anchor_bin=ANCHOR, n_mc=n_mc, seed=seed,
        f_lo=F_LO, f_hi=F_HI, xtol=XTOL, include_wobble=True, n_bootstrap_final=0, return_raw=True)


# ----------------------------------------------------------------------------------------
# Estatisticas
# ----------------------------------------------------------------------------------------

def _bin_deltas(lr_real, lr_mock, bin_idx, in_range, res, mc):
    rp = lr_real[mc, res]
    mp = lr_mock[mc, res]
    bi = bin_idx[res]
    ir = in_range[res]
    out = np.full(NB, np.nan)
    for b in range(NB):
        m = ir & (bi == b)
        if m.any():
            out[b] = np.median(rp[m]) - np.median(mp[m])
    return out


def fit_a0_joint(gN: np.ndarray, delta: np.ndarray, w: np.ndarray, x0: float = 1.0):
    """Ajuste de a0 unico a (K x NB) pontos, pesos = contagens (peso deterministico declarado)."""
    gN = np.asarray(gN, float).ravel()
    d = np.asarray(delta, float).ravel()
    w = np.asarray(w, float).ravel()
    ok = np.isfinite(d) & np.isfinite(gN) & (w > 0)
    if ok.sum() < 3:
        return None
    sigma = 1.0 / np.sqrt(w[ok] / w[ok].mean())

    def model(g, x):
        return sp.delta_aqual(g, x * 1e-10)
    try:
        popt, _ = curve_fit(model, gN[ok], d[ok], p0=[x0], sigma=sigma, maxfev=20000)
        a0 = float(popt[0]) * 1e-10
        if not np.isfinite(a0) or a0 <= 0:
            return None
        return a0
    except Exception:
        return None


def verdict(ci_lo, ci_hi):
    if ci_lo is None or ci_hi is None:
        return "NO_FIT"
    a_in = ci_lo <= A0_A <= ci_hi
    b_in = ci_lo <= A0_B <= ci_hi
    if a_in and b_in:
        return "BOTH_SURVIVE"
    if a_in:
        return "A_ONLY"
    if b_in:
        return "B_ONLY"
    return "BOTH_FALSIFIED"


def pct_ci(arr, q=(2.5, 97.5)):
    a = np.asarray([x for x in arr if x is not None and np.isfinite(x)], float)
    if len(a) < 20:
        return None, None
    lo, hi = np.percentile(a, q)
    return float(lo), float(hi)


def fpb(strata: list, n_boot: int, rng: np.random.Generator, recal: bool):
    """Bootstrap do PROCEDIMENTO completo por estrato (reamostra sistemas e realizacoes MC;
    recalibra f por interpolacao/extrapolacao LINEAR em f entre dois cortes pareados
    (raw_a em f_a=f_hat, raw_b em f_b=f_a+-H_F, mesma seed => sorteios nested).
    Retorna f_rep (n_boot,K) e delta (n_boot,K,NB)."""
    K = len(strata)
    f_rep = np.full((n_boot, K), np.nan)
    delta = np.full((n_boot, K, NB), np.nan)
    for rep in range(n_boot):
        for k, st in enumerate(strata):
            n = st["n_sys"]
            res = rng.integers(0, n, size=n)
            mc = rng.integers(0, st["n_mc"], size=n)
            da = _bin_deltas(st["raw_a"]["log_ratio_real"], st["raw_a"]["log_ratio_mock"],
                             st["raw_a"]["bin_idx"], st["raw_a"]["in_range"], res, mc)
            if not recal:
                delta[rep, k] = da
                f_rep[rep, k] = st["f_a"]
                continue
            db = _bin_deltas(st["raw_b"]["log_ratio_real"], st["raw_b"]["log_ratio_mock"],
                             st["raw_b"]["bin_idx"], st["raw_b"]["in_range"], res, mc)
            dfb = st["f_b"] - st["f_a"]
            slope = (db[ANCHOR] - da[ANCHOR]) / dfb
            if not np.isfinite(slope) or abs(slope) < 1e-9:
                continue
            f_new = st["f_a"] - da[ANCHOR] / slope
            f_new = min(max(f_new, F_LO), F_HI)
            t = (f_new - st["f_a"]) / dfb
            delta[rep, k] = da + t * (db - da)
            f_rep[rep, k] = f_new
    return f_rep, delta


def anchor_resid_bootstrap(raw: dict, label: np.ndarray, K: int, n_boot: int, rng: np.random.Generator):
    """Resido do bin-ancora por estrato SOB f escalar fixo (teste D de heterogeneidade)."""
    n = raw["n_sys"]
    n_mc = raw["n_mc"]
    lr, lm = raw["log_ratio_real"], raw["log_ratio_mock"]
    bi, ir = raw["bin_idx"], raw["in_range"]
    point = []
    for k in range(K):
        mk = ir & (bi == ANCHOR) & (label == k)
        point.append(float(np.median(lr[:, mk]) - np.median(lm[:, mk])))
    boots = np.full((n_boot, K), np.nan)
    for rep in range(n_boot):
        res = rng.integers(0, n, size=n)
        mc = rng.integers(0, n_mc, size=n)
        rp, mp = lr[mc, res], lm[mc, res]
        bir, irr, lab = bi[res], ir[res], label[res]
        for k in range(K):
            m = irr & (bir == ANCHOR) & (lab == k)
            if m.any():
                boots[rep, k] = np.median(rp[m]) - np.median(mp[m])
    return np.array(point), boots


# ----------------------------------------------------------------------------------------
# Uma replicata completa
# ----------------------------------------------------------------------------------------

def run_replicate(task: dict) -> dict:
    t0 = time.time()
    cfg_id, rep, cfg = task["config_id"], task["rep"], task["cfg"]
    n_mc, n_boot = task.get("n_mc", N_MC), task.get("n_boot", N_BOOT)
    seeds = task.get("seeds") or {r: seed_for(cfg_id, rep, i) for i, r in enumerate(
        ["cat", "truth", "calib_scalar", "calib_low", "calib_high", "boot"])}

    n_sys = cfg["n_sys"]
    cat = generate_equal_bin_catalog(n_sys, np.random.default_rng(seeds["cat"]))
    n_sys = cat["n_sys"]
    rng_truth = np.random.default_rng(seeds["truth"])
    if cfg.get("pi_bin") is None:
        pi_vec = np.full(n_sys, PI_HIGH)
    else:   # fracao RUWE-alto depende do bin de gN (cenarios hipoteticos H2)
        lg_cat = np.log10(sp.G_SI * cat["Mtot_cat_kg"] / cat["s_m"] ** 2)
        pi_vec = np.asarray(cfg["pi_bin"])[np.digitize(lg_cat, BIN_EDGES[1:-1])]
    label = (rng_truth.random(n_sys) < pi_vec).astype(int)      # 0 = RUWE-baixo, 1 = RUWE-alto
    f_true_vec = np.where(label == 1, cfg["f_high"], cfg["f_low"])
    wob = np.where(label == 1, cfg["wob"], 1.0)
    kap = np.where(label == 1, cfg["kap"], 1.0)
    ext = np.where(label == 1, cfg.get("ext_kms", 0.0) * 1000.0, 0.0) if cfg.get("ext_kms", 0.0) else None
    truth = build_truth(cat, f_true_vec, cfg["a0"], rng_truth, wobble_scale=wob, noise_kappa=kap, extra_v_si=ext)
    v_p = truth["v_p"]
    a0_true = cfg["a0"]

    out = dict(config_id=cfg_id, rep=rep, name=cfg["name"], n_sys=n_sys, seeds=seeds,
               frac_high=float(label.mean()), frac_has_multi_true=float(truth["has_multi"].mean()))

    # ---------------- Metodo ESCALAR (Estagio 2 travado) ----------------
    all_idx = np.arange(n_sys)
    cs = calibrate(cat, v_p, seeds["calib_scalar"], n_mc)
    f_s = cs["f_multi_calibrated"]
    fr = cs["final_result"]
    delta_s = np.array(fr["delta_obs_newt_primary"])
    gN_s = np.array(fr["gN_bin_median_si"])
    a0_s = sp.fit_a0(gN_s, delta_s, x0=1.0)
    ci_s = sp.bootstrap_a0_refit(cs["final_raw"], NB, n_boot, seeds["boot"], x0_list=(1.0,))["x0=1.0"]
    out["scalar"] = dict(f_hat=f_s, converged=cs["converged_bracket"], delta=delta_s.tolist(),
                         a0=a0_s, ci_lo=ci_s["ci95_lo_si_m_s2"], ci_hi=ci_s["ci95_hi_si_m_s2"],
                         verdict=verdict(ci_s["ci95_lo_si_m_s2"], ci_s["ci95_hi_si_m_s2"]),
                         n_per_bin=fr["n_systems_per_bin"])

    # ---------------- Teste D: heterogeneidade do residuo do ancora por estrato ----------------
    rng_b = np.random.default_rng(seeds["boot"])
    point_D, boots_D = anchor_resid_bootstrap(cs["final_raw"], label, 2, n_boot, rng_b)
    Dpt = float(point_D[1] - point_D[0])
    Dboot = boots_D[:, 1] - boots_D[:, 0]
    D_lo, D_hi = pct_ci(Dboot)
    out["D_test"] = dict(delta_anchor_low=float(point_D[0]), delta_anchor_high=float(point_D[1]),
                         D=Dpt, D_ci_lo=D_lo, D_ci_hi=D_hi,
                         detects_hetero_95=bool(D_lo is not None and D_lo > 0),
                         sd_D=float(np.nanstd(Dboot)))

    # ---------------- Metodo CONDICIONAL: calibracao por estrato ----------------
    strata = []
    cond = {}
    for k, role in ((0, "calib_low"), (1, "calib_high")):
        idx = np.where(label == k)[0]
        cs_k = cat_k = None
        cat_k = subset_cat(cat, idx)
        calib_k = calibrate(cat_k, v_p[idx], seeds[role], n_mc)
        f_k = calib_k["f_multi_calibrated"]
        f_b = f_k + H_F if f_k + H_F <= F_HI else f_k - H_F
        raw_b = run_pipeline(cat_k, v_p[idx], f_b, seeds[role] + 777, n_mc)
        frk = calib_k["final_result"]
        strata.append(dict(n_sys=cat_k["n_sys"], n_mc=n_mc, raw_a=calib_k["final_raw"], raw_b=raw_b["_raw"],
                           f_a=f_k, f_b=f_b))
        cond[k] = dict(f_hat=f_k, converged=calib_k["converged_bracket"],
                       boundary_zero=bool((not calib_k["converged_bracket"]) and f_k <= F_LO + 1e-12),
                       boundary_top=bool((not calib_k["converged_bracket"]) and f_k >= F_HI - 1e-12),
                       delta=frk["delta_obs_newt_primary"], gN=frk["gN_bin_median_si"],
                       n_per_bin=frk["n_systems_per_bin"], n_sys=cat_k["n_sys"],
                       n_anchor=frk["n_systems_per_bin"][ANCHOR])
        a0_k = sp.fit_a0(np.array(frk["gN_bin_median_si"]), np.array(frk["delta_obs_newt_primary"]), x0=1.0)
        cond[k]["a0_stratum"] = a0_k
    out["cond_strata"] = cond
    f_hat_vec = np.where(label == 1, cond[1]["f_hat"], cond[0]["f_hat"])
    out["f_true"] = dict(low=cfg["f_low"], high=cfg["f_high"])

    # estimativa pontual conjunta (ajuste unico de a0 a K x NB pontos, pesos = contagens)
    gN_k = np.array([cond[k]["gN"] for k in (0, 1)])
    d_k = np.array([cond[k]["delta"] for k in (0, 1)])
    w_k = np.array([cond[k]["n_per_bin"] for k in (0, 1)], float)
    a0_c = fit_a0_joint(gN_k, d_k, w_k)

    # bootstrap do procedimento (com e sem recalibracao de f)
    res_fpb = {}
    for tag, recal in (("recal", True), ("fixed", False)):
        f_rep, delta_rep = fpb(strata, n_boot, np.random.default_rng(seeds["boot"] + (1 if recal else 2)), recal)
        a0_reps = []
        for rep_i in range(n_boot):
            if np.any(~np.isfinite(delta_rep[rep_i])):
                continue
            a0_reps.append(fit_a0_joint(gN_k, delta_rep[rep_i], w_k))
        lo, hi = pct_ci(a0_reps)
        sd_kb = np.nanstd(delta_rep, axis=0)  # (K, NB)
        res_fpb[tag] = dict(ci_lo=lo, ci_hi=hi, verdict=verdict(lo, hi),
                            n_valid=int(sum(x is not None for x in a0_reps)),
                            sd_kb=sd_kb.tolist(),
                            f_rep_sd=np.nanstd(f_rep, axis=0).tolist(),
                            f_rep_ci=[list(pct_ci(f_rep[:, k])) for k in (0, 1)])
        if recal:
            # diferenca de f entre estratos, com IC do bootstrap do procedimento (estatistica que o
            # metodo CONDICIONAL usa para distinguir heterogeneo de homogeneo)
            df_rep = f_rep[:, 1] - f_rep[:, 0]
            df_lo, df_hi = pct_ci(df_rep)
            res_fpb[tag]["df_point"] = float(cond[1]["f_hat"] - cond[0]["f_hat"])
            res_fpb[tag]["df_ci"] = [df_lo, df_hi]
            res_fpb[tag]["df_detects_hetero_95"] = bool(df_lo is not None and df_lo > 0)
            Delta_rep = delta_rep[:, 1, :4] - delta_rep[:, 0, :4]   # bins 0..3 (ancora excluido)
            res_fpb[tag]["Delta_sd"] = np.nanstd(Delta_rep, axis=0).tolist()
            Dpoint = (d_k[1, :4] - d_k[0, :4])
            sdD = np.nanstd(Delta_rep, axis=0)
            z = Dpoint / np.where(sdD > 0, sdD, np.nan)
            T = float(np.nansum(z ** 2))
            res_fpb[tag]["Delta_point"] = Dpoint.tolist()
            res_fpb[tag]["z_bins0_3"] = z.tolist()
            res_fpb[tag]["T_univ"] = T
            res_fpb[tag]["T_univ_flag_chi2_4_95"] = bool(T > stats.chi2.ppf(0.95, 4))
            res_fpb[tag]["zmax_flag_bonf4"] = bool(np.nanmax(np.abs(z)) > stats.norm.ppf(1 - 0.05 / 8))
            if a0_c is not None:
                model = sp.delta_aqual(gN_k[:, :4], a0_c)
                r = (d_k[:, :4] - model) / np.where(sd_kb[:, :4] > 0, sd_kb[:, :4], np.nan)
                chi2 = float(np.nansum(r ** 2))
                dof = 2 * 4 - 1
                res_fpb[tag]["gof_chi2"] = chi2
                res_fpb[tag]["gof_dof"] = dof
                res_fpb[tag]["gof_flag_95"] = bool(chi2 > stats.chi2.ppf(0.95, dof))
    out["cond_joint"] = dict(a0=a0_c, **{k: v for k, v in res_fpb.items()})

    # alternativa "corte RUWE<1,4": usa so' o estrato baixo com o metodo escalar travado
    st0 = strata[0]
    ci_low = sp.bootstrap_a0_refit(st0["raw_a"], NB, n_boot, seeds["boot"] + 3, x0_list=(1.0,))["x0=1.0"]
    out["cut_low_only"] = dict(a0=cond[0]["a0_stratum"], ci_lo=ci_low["ci95_lo_si_m_s2"],
                               ci_hi=ci_low["ci95_hi_si_m_s2"],
                               verdict=verdict(ci_low["ci95_lo_si_m_s2"], ci_low["ci95_hi_si_m_s2"]),
                               f_hat=cond[0]["f_hat"], n_sys=cond[0]["n_sys"])

    # aproximacao linear do bootstrap: erro vs. rodada exata em f_rep da amostra original
    out["a0_true"] = a0_true
    out["elapsed_s"] = time.time() - t0
    return out


# ----------------------------------------------------------------------------------------
# Testes unitarios (Parte "unit")
# ----------------------------------------------------------------------------------------

def part_unit(log) -> dict:
    res = {}
    n = 1500
    cat = bp.generate_structural_catalog(n, np.random.default_rng(dev_seed(0)))
    # U1: equivalencia exata build_truth(f escalar) vs bp.build_fake_real_dataset
    ref = bp.build_fake_real_dataset(cat, 0.3, A0_A, True, np.random.default_rng(dev_seed(1)))
    mine = build_truth(cat, np.full(n, 0.3), A0_A, np.random.default_rng(dev_seed(1)))
    d = float(np.max(np.abs(ref["v_p_fake_real"] - mine["v_p"])))
    res["U1_truth_equivalence_max_abs_diff"] = d
    res["U1_pass"] = bool(d == 0.0)
    log(f"U1 build_truth(f escalar) == bp.build_fake_real_dataset: max|dv|={d:.3e} -> {'PASS' if d == 0.0 else 'FAIL'} (seed {dev_seed(1)})")

    # U2: pipeline com f vetorial constante == pipeline escalar
    v_p = mine["v_p"]
    r_sc = run_pipeline(cat, v_p, 0.2, dev_seed(2), 40)
    r_vec = run_pipeline(cat, v_p, 0.2, dev_seed(2), 40, f_vec=np.full(n, 0.2))
    d2 = float(np.nanmax(np.abs(np.array(r_sc["delta_obs_newt_primary"]) - np.array(r_vec["delta_obs_newt_primary"]))))
    res["U2_fvector_equals_scalar_max_abs_diff"] = d2
    res["U2_pass"] = bool(d2 == 0.0)
    log(f"U2 pipeline f vetorial constante == escalar: max|d(delta)|={d2:.3e} -> {'PASS' if d2 == 0.0 else 'FAIL'} (seed {dev_seed(2)})")

    # U3: aproximacao linear de f no bootstrap: delta interpolado linearmente vs rodada exata
    f_a, f_b, f_mid = 0.10, 0.25, 0.17
    ra = run_pipeline(cat, v_p, f_a, dev_seed(3), 60)
    rb = run_pipeline(cat, v_p, f_b, dev_seed(3), 60)
    rm = run_pipeline(cat, v_p, f_mid, dev_seed(3), 60)
    da, db, dm = (np.array(x["delta_obs_newt_primary"]) for x in (ra, rb, rm))
    t = (f_mid - f_a) / (f_b - f_a)
    lin = da + t * (db - da)
    res["U3_linear_interp_err_per_bin"] = (lin - dm).tolist()
    res["U3_max_abs_err"] = float(np.nanmax(np.abs(lin - dm)))
    log(f"U3 erro da interpolacao linear em f (f=0,10..0,25, ponto 0,17; n={n}, n_mc=60; seed {dev_seed(3)}): "
        f"max={res['U3_max_abs_err']:.4f} dex; por bin={np.round(lin - dm, 4).tolist()}")

    # U4: sanidade do pipeline: populacao homogenea f=0,30, recupera f escalar ~0,30
    cat2 = bp.generate_structural_catalog(6000, np.random.default_rng(dev_seed(4)))
    tr2 = build_truth(cat2, np.full(6000, 0.30), None, np.random.default_rng(dev_seed(5)))
    c2 = calibrate(cat2, tr2["v_p"], dev_seed(6), 60)
    res["U4_homog_f0.30_recovered"] = c2["f_multi_calibrated"]
    res["U4_converged"] = c2["converged_bracket"]
    log(f"U4 populacao homogenea f_true=0,30, N=6000, n_mc=60: f_hat={c2['f_multi_calibrated']:.4f} "
        f"converged={c2['converged_bracket']} (seeds {dev_seed(4)},{dev_seed(5)},{dev_seed(6)})")
    nb = [int(((np.log10(sp.G_SI * cat2['Mtot_cat_kg'] / cat2['s_m'] ** 2) >= BIN_EDGES[b]) &
               (np.log10(sp.G_SI * cat2['Mtot_cat_kg'] / cat2['s_m'] ** 2) < BIN_EDGES[b + 1])).sum()) for b in range(NB)]
    res["catalog_counts_per_bin_N6000"] = nb
    log(f"   contagens por bin do catalogo sintetico (N=6000): {nb}  (o real tem ~6040 por bin com N=30203)")
    # U5: catalogo de contagens iguais por bin
    cat3 = generate_equal_bin_catalog(6000, np.random.default_rng(dev_seed(7)))
    lg3 = np.log10(sp.G_SI * cat3["Mtot_cat_kg"] / cat3["s_m"] ** 2)
    nb3 = np.bincount(np.digitize(lg3, BIN_EDGES[1:-1]), minlength=NB).tolist()
    res["U5_equal_bin_counts_N6000"] = nb3
    log(f"U5 catalogo de contagens iguais por bin (N=6000, seed {dev_seed(7)}): {nb3}")
    return res


def _explore_work(args):
    label, wob, kap, ext = args
    seeds = {r: dev_seed(8 + i) for i, r in enumerate(
        ["cat", "truth", "calib_scalar", "calib_low", "calib_high", "boot"])}   # 968..973, reusadas em todos os botoes
    cfg = dict(name=label, f_low=0.02, f_high=0.45, a0=A0_A, n_sys=N_HOLDOUT_DOCUMENTED, wob=wob, kap=kap,
               ext_kms=ext, pi_bin=None)
    try:
        r = run_replicate(dict(config_id=0, rep=0, cfg=cfg, n_mc=N_MC, n_boot=200, seeds=seeds))
    except Exception as e:  # noqa
        import traceback
        return dict(label=label, error=repr(e), tb=traceback.format_exc())
    rc = r["cond_joint"]["recal"]
    return dict(label=label, wobble_scale_high=wob, noise_kappa_high=kap, ext_kms_high=ext,
                f_hat_low=r["cond_strata"][0]["f_hat"], f_hat_high=r["cond_strata"][1]["f_hat"],
                D_scalar_anchor=r["D_test"]["D"], Delta_bins0_3=rc["Delta_point"], z_bins0_3=rc["z_bins0_3"],
                T_univ=rc["T_univ"], gof_chi2=rc.get("gof_chi2"),
                a0_scalar=r["scalar"]["a0"], a0_cond_joint=r["cond_joint"]["a0"], a0_cut_low=r["cut_low_only"]["a0"],
                a0_stratum_low=r["cond_strata"][0]["a0_stratum"], a0_stratum_high=r["cond_strata"][1]["a0_stratum"],
                delta_low=r["cond_strata"][0]["delta"], delta_high=r["cond_strata"][1]["delta"], seeds=seeds)


def part_explore(log, workers: int) -> dict:
    """Exploracao de BOTOES de mis-especificacao do estrato RUWE-alto (1 replicata, N=12944, seeds de
    dev 968-973 reusadas em todos os botoes => dados-base identicos, so' o botao muda). Serve para
    escolher o valor do artefato usado nas configs 3, 4 e 13 e documentar o que NAO funcionou."""
    import multiprocessing as mp
    combos = [("wobble_x4", 4.0, 1.0, 0.0), ("wobble_x8", 8.0, 1.0, 0.0), ("wobble_x16", 16.0, 1.0, 0.0),
              ("noise_x3", 1.0, 3.0, 0.0), ("noise_x5", 1.0, 5.0, 0.0), ("noise_x8", 1.0, 8.0, 0.0)] + \
             [(f"artifact_{e}kms", 1.0, 1.0, e) for e in (0.05, 0.1, 0.2, 0.3, 0.5, 0.8, 1.2, 2.0)]
    with mp.get_context("fork").Pool(workers) as pool:
        outs = pool.map(_explore_work, combos)
    for o in outs:
        if "error" in o:
            log(f"EXPLORE ERRO {o['label']}: {o['error']}")
            continue
        log(f"EXPLORE {o['label']:>16}: f_hat(low,high)=({o['f_hat_low']:.3f},{o['f_hat_high']:.3f}) D={o['D_scalar_anchor']:+.3f} "
            f"Delta_b0-3={np.round(o['Delta_bins0_3'], 3).tolist()} T={o['T_univ']:.1f} | a0 escalar={o['a0_scalar']:.3g} "
            f"cond={o['a0_cond_joint']:.3g} cut_low={o['a0_cut_low']:.3g} (a0_true={A0_A:.3g})")
    return {"combos": outs}


# ----------------------------------------------------------------------------------------
# Retrospectiva aritmetica: GOF da forma AQUAL sobre os numeros JA PUBLICADOS do Estagio 2
# (nenhum dado real; so' aritmetica sobre valores citados)
# ----------------------------------------------------------------------------------------

def part_retro_gof(log) -> dict:
    # RESULTS_PRIMARY_STAGE2.md:95-99 (delta corrigido e IC95% por bin, N_bootstrap=2000)
    delta = np.array([0.0777, 0.0657, 0.0486, 0.0417, -0.0002])
    lo = np.array([0.0495, 0.0399, 0.0223, 0.0164, -0.0254])
    hi = np.array([0.1041, 0.0924, 0.0746, 0.0686, 0.0255])
    gN = np.array([1.8807e-10, 1.6645e-9, 6.2203e-9, 1.7190e-8, 4.8638e-8])  # RESULTS_PRIMARY_STAGE2.md:95-99
    sigma = (hi - lo) / (2 * 1.96)   # aproximacao gaussiana do IC95% (ressalva: ignora covariancias e incerteza de f_hat)
    out = {"sigma_approx": sigma.tolist()}
    for name, a0 in (("a0_fit_stage2=6.125e-11 (RESULTS_PRIMARY_STAGE2.md:32)", 6.125e-11),
                     ("a0_A=1.082288e-10", A0_A), ("a0_B=6.800218e-10", A0_B)):
        model = sp.delta_aqual(gN, a0)
        r = (delta - model) / sigma
        chi2_all = float(np.sum(r ** 2))
        chi2_nonanchor = float(np.sum(r[:4] ** 2))
        out[name] = dict(model=np.round(model, 5).tolist(), resid_over_sigma=np.round(r, 2).tolist(),
                         chi2_5bins=chi2_all, chi2_bins0_3=chi2_nonanchor,
                         p_5bins_dof4=float(stats.chi2.sf(chi2_all, 4)))
        log(f"RETRO-GOF {name}: modelo={np.round(model, 4).tolist()} chi2(5 bins, dof~4)={chi2_all:.1f} "
            f"p={stats.chi2.sf(chi2_all, 4):.2e} | chi2(bins 0-3)={chi2_nonanchor:.1f}")
    log("   ressalva: sigma = (hi-lo)/3.92 do IC95% publicado; ignora covariancias e incerteza de f_hat; "
        "serve so' como ordem de grandeza.")

    # --- universalidade RETROSPECTIVA entre estratos RUWE, com a tabela do debunker ---
    # ADVERSARIAL_DEBUNKER_REPORT.md:319-323 (corrigido com f_hat escalar=0,1037, N_bootstrap=500)
    hi_d = np.array([0.486, 0.458, 0.412, 0.325, 0.171])
    hi_lo = np.array([0.404, 0.392, 0.348, 0.265, 0.126]); hi_hi = np.array([0.559, 0.531, 0.476, 0.391, 0.212])
    lo_d = np.array([0.020, 0.002, -0.012, -0.014, -0.066])
    lo_lo = np.array([-0.008, -0.026, -0.042, -0.042, -0.096]); lo_hi = np.array([0.050, 0.033, 0.015, 0.013, -0.038])
    s_h = (hi_hi - hi_lo) / 3.92
    s_l = (lo_hi - lo_lo) / 3.92
    # ancorando cada estrato no proprio bin 4 (o que a calibracao condicional faria; aproximacao: perfil ~ invariante em f)
    Db = (hi_d[:4] - hi_d[4]) - (lo_d[:4] - lo_d[4])
    sD = np.sqrt(s_h[:4] ** 2 + s_h[4] ** 2 + s_l[:4] ** 2 + s_l[4] ** 2)
    z = Db / sD
    T = float(np.sum(z ** 2))
    out["retro_universality_debunker_table"] = dict(
        anchored_high=(hi_d[:4] - hi_d[4]).round(3).tolist(), anchored_low=(lo_d[:4] - lo_d[4]).round(3).tolist(),
        Delta_bins0_3=Db.round(3).tolist(), sigma=sD.round(3).tolist(), z=z.round(1).tolist(), T=T,
        p_chi2_4=float(stats.chi2.sf(T, 4)))
    log(f"RETRO-UNIVERSALIDADE (tabela do debunker, ADVERSARIAL_DEBUNKER_REPORT.md:319-323): "
        f"perfil ancorado alto={ (hi_d[:4]-hi_d[4]).round(3).tolist() } baixo={ (lo_d[:4]-lo_d[4]).round(3).tolist() } "
        f"Delta_b={Db.round(3).tolist()} z={z.round(1).tolist()} T={T:.1f} (p chi2_4={stats.chi2.sf(T,4):.1e})")
    log("   ressalva: aritmetica gaussiana sobre IC95% (N_boot=500) publicados; assume perfil ancorado ~ invariante a f; "
        "NAO e' um teste pre-registrado nem uma analise nova -- o debunker ja' viu a amostra de descoberta.")
    return out


# ----------------------------------------------------------------------------------------
# Agregacao
# ----------------------------------------------------------------------------------------

def wilson(k, n, z=1.96):
    if n == 0:
        return [None, None]
    p = k / n
    den = 1 + z * z / n
    c = (p + z * z / (2 * n)) / den
    h = z * math.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / den
    return [round(c - h, 3), round(c + h, 3)]


def _bias(a0, a0_true):
    if a0 is None or a0_true is None or a0 <= 0:
        return None
    return math.log10(a0 / a0_true)


def summarize(config_id: int, results: list) -> dict:
    cfg = CONFIGS[config_id]
    a0_true = cfg["a0"]
    n = len(results)
    S = dict(config_id=config_id, name=cfg["name"], n_reps=n, n_sys=cfg["n_sys"],
             f_true=dict(low=cfg["f_low"], high=cfg["f_high"]), a0_true=a0_true,
             wobble_scale_high=cfg["wob"], noise_kappa_high=cfg["kap"], artifact_kms_high=cfg.get("ext_kms", 0.0))

    def stat(vals):
        v = np.array([x for x in vals if x is not None and np.isfinite(x)], float)
        if len(v) == 0:
            return dict(n=0)
        return dict(n=int(len(v)), mean=float(v.mean()), sd=float(v.std(ddof=1)) if len(v) > 1 else None,
                    median=float(np.median(v)), min=float(v.min()), max=float(v.max()))

    S["n_anchor_low"] = stat([r["cond_strata"][0]["n_anchor"] for r in results])
    S["n_anchor_high"] = stat([r["cond_strata"][1]["n_anchor"] for r in results])
    S["n_sys_high"] = stat([r["cond_strata"][1]["n_sys"] for r in results])
    S["n_bin0_high"] = stat([r["cond_strata"][1]["n_per_bin"][0] for r in results])
    S["f_hat_scalar"] = stat([r["scalar"]["f_hat"] for r in results])
    S["f_hat_low"] = stat([r["cond_strata"][0]["f_hat"] for r in results])
    S["f_hat_high"] = stat([r["cond_strata"][1]["f_hat"] for r in results])
    S["boundary_zero_rate_low"] = float(np.mean([r["cond_strata"][0]["boundary_zero"] for r in results]))
    S["boundary_zero_rate_high"] = float(np.mean([r["cond_strata"][1]["boundary_zero"] for r in results]))
    S["not_converged_rate_scalar"] = float(np.mean([not r["scalar"]["converged"] for r in results]))
    S["df_point"] = stat([r["cond_joint"]["recal"]["df_point"] for r in results])
    kk = sum(bool(r["cond_joint"]["recal"]["df_detects_hetero_95"]) for r in results)
    S["df_detection_rate_95"] = dict(rate=kk / n, k=kk, n=n, wilson95=wilson(kk, n))
    S["D_stat"] = stat([r["D_test"]["D"] for r in results])
    S["D_sd_boot"] = stat([r["D_test"]["sd_D"] for r in results])
    k = sum(r["D_test"]["detects_hetero_95"] for r in results)
    S["D_detection_rate_95"] = dict(rate=k / n, k=k, n=n, wilson95=wilson(k, n))

    methods = {
        "scalar": lambda r: (r["scalar"]["a0"], r["scalar"]["ci_lo"], r["scalar"]["ci_hi"], r["scalar"]["verdict"]),
        "cond_joint_recal": lambda r: (r["cond_joint"]["a0"], r["cond_joint"]["recal"]["ci_lo"], r["cond_joint"]["recal"]["ci_hi"], r["cond_joint"]["recal"]["verdict"]),
        "cond_joint_fixedf": lambda r: (r["cond_joint"]["a0"], r["cond_joint"]["fixed"]["ci_lo"], r["cond_joint"]["fixed"]["ci_hi"], r["cond_joint"]["fixed"]["verdict"]),
        "cut_low_only": lambda r: (r["cut_low_only"]["a0"], r["cut_low_only"]["ci_lo"], r["cut_low_only"]["ci_hi"], r["cut_low_only"]["verdict"]),
    }
    S["methods"] = {}
    for mname, fn in methods.items():
        rows = [fn(r) for r in results]
        biases = [_bias(a, a0_true) for a, _, _, _ in rows]
        m = dict(bias_log10=stat(biases),
                 n_no_fit=int(sum(a is None for a, _, _, _ in rows)))
        ci_w = [math.log10(hi / lo) if (lo and hi and lo > 0) else None for _, lo, hi, _ in rows]
        m["ci_width_log10"] = stat(ci_w)
        if a0_true is not None:
            cov = [(lo is not None and lo <= a0_true <= hi) for _, lo, hi, _ in rows if lo is not None]
            m["coverage_a0_true"] = dict(rate=(sum(cov) / len(cov)) if cov else None, k=int(sum(cov)), n=len(cov),
                                         wilson95=wilson(sum(cov), len(cov)))
        vc = {}
        for _, _, _, v in rows:
            vc[v] = vc.get(v, 0) + 1
        m["verdict_counts"] = vc
        S["methods"][mname] = m

    # gates (somente metodo condicional, recal)
    def rate(key_fn):
        vals = [key_fn(r) for r in results]
        vals = [v for v in vals if v is not None]
        return dict(rate=(sum(vals) / len(vals)) if vals else None, k=int(sum(vals)), n=len(vals),
                    wilson95=wilson(int(sum(vals)), len(vals)))
    S["gate_universality_T_chi2_4_95"] = rate(lambda r: r["cond_joint"]["recal"].get("T_univ_flag_chi2_4_95"))
    S["gate_universality_zmax_bonf"] = rate(lambda r: r["cond_joint"]["recal"].get("zmax_flag_bonf4"))
    S["gate_gof_95"] = rate(lambda r: r["cond_joint"]["recal"].get("gof_flag_95"))
    S["T_univ"] = stat([r["cond_joint"]["recal"].get("T_univ") for r in results])
    S["gof_chi2"] = stat([r["cond_joint"]["recal"].get("gof_chi2") for r in results])
    S["delta_high_minus_low_bins0_3_mean"] = np.nanmean(
        [r["cond_joint"]["recal"]["Delta_point"] for r in results], axis=0).tolist()
    S["f_rep_sd_high_mean"] = float(np.nanmean([r["cond_joint"]["recal"]["f_rep_sd"][1] for r in results]))
    S["f_rep_sd_low_mean"] = float(np.nanmean([r["cond_joint"]["recal"]["f_rep_sd"][0] for r in results]))
    S["elapsed_s_mean"] = float(np.mean([r["elapsed_s"] for r in results]))
    S["gN_bin_median_synthetic_low_stratum_mean"] = np.nanmean(
        [r["cond_strata"][0]["gN"] for r in results], axis=0).tolist()
    S["stratum_a0_log10_ratio_high_over_low"] = stat([
        (math.log10(r["cond_strata"][1]["a0_stratum"] / r["cond_strata"][0]["a0_stratum"])
         if (r["cond_strata"][1]["a0_stratum"] and r["cond_strata"][0]["a0_stratum"]) else None) for r in results])
    S["stratum_a0_high"] = stat([r["cond_strata"][1]["a0_stratum"] for r in results])
    S["stratum_a0_low"] = stat([r["cond_strata"][0]["a0_stratum"] for r in results])
    # criterios de identificabilidade (limiares DECLARADOS aqui, antes de ver o resultado final)
    LOG_AB = math.log10(A0_B / A0_A)
    crit = dict(C1_detect_hetero_power_ge_0p80_or_fp_le_0p10=None,
                C2_cond_abs_median_bias_le_0p30dex=None,
                C3_cond_ci_width_median_lt_log10_aB_over_aA=None, log10_aB_over_aA=LOG_AB)
    if cfg["f_low"] != cfg["f_high"]:
        crit["C1_detect_hetero_power_ge_0p80_or_fp_le_0p10"] = bool(
            S["D_detection_rate_95"]["rate"] >= 0.80 and S["df_detection_rate_95"]["rate"] >= 0.80)
    else:
        crit["C1_detect_hetero_power_ge_0p80_or_fp_le_0p10"] = bool(
            S["D_detection_rate_95"]["rate"] <= 0.10 and S["df_detection_rate_95"]["rate"] <= 0.10)
    mb = S["methods"]["cond_joint_recal"]["bias_log10"]
    if mb.get("n", 0) > 0:
        crit["C2_cond_abs_median_bias_le_0p30dex"] = bool(abs(mb["median"]) <= 0.30)
    mw = S["methods"]["cond_joint_recal"]["ci_width_log10"]
    if mw.get("n", 0) > 0:
        crit["C3_cond_ci_width_median_lt_log10_aB_over_aA"] = bool(mw["median"] < LOG_AB)
    S["criteria"] = crit
    return S


STAGE2_BIAS_F_LO, STAGE2_BIAS_F_HI = 1.4, 1.6  # PREREGISTRATION_STAGE2.md:341-343, RESULTS_PRIMARY_STAGE2.md:159-163


def decide_stage2(r: dict) -> str:
    """Regra do Estagio 2 (travada): mecanica (IC do escalar) + camada interpretativa (CI/1,6 .. CI/1,4)."""
    lo, hi = r["scalar"]["ci_lo"], r["scalar"]["ci_hi"]
    if lo is None:
        return "NO_FIT"
    v = verdict(lo, hi)
    if v == "BOTH_SURVIVE":
        return "NAO_DISTINGUE_CI"
    v2 = verdict(lo / STAGE2_BIAS_F_HI, hi / STAGE2_BIAS_F_LO)
    return v if v2 == v else "NAO_DISTINGUE_VIES"


def decide_stage3(r: dict, thr: dict, F3: float) -> str:
    """Regra PROPOSTA (rascunho): gates (universalidade, GOF) -> CI conjunto -> camada de vies bilateral."""
    rc = r["cond_joint"]["recal"]
    lo, hi = rc["ci_lo"], rc["ci_hi"]
    if lo is None:
        return "NO_FIT"
    zmax = float(np.max(np.abs(np.array(rc["z_bins0_3"]))))
    gate_fail = (rc["T_univ"] > thr["T_univ"]) or (zmax > thr["zmax"]) or (rc.get("gof_chi2", 0.0) > thr["gof_chi2"])
    if gate_fail:
        return "NAO_DISTINGUE_GATE"
    v = verdict(lo, hi)
    if v == "BOTH_SURVIVE":
        return "NAO_DISTINGUE_CI"
    v2 = verdict(lo / F3, hi * F3)
    return v if v2 == v else "NAO_DISTINGUE_VIES"


def decide_cut_low(r: dict, F3: float) -> str:
    """Alternativa 'corte RUWE<=1,4': so' o estrato baixo, metodo escalar travado, camada de vies bilateral."""
    lo, hi = r["cut_low_only"]["ci_lo"], r["cut_low_only"]["ci_hi"]
    if lo is None:
        return "NO_FIT"
    v = verdict(lo, hi)
    if v == "BOTH_SURVIVE":
        return "NAO_DISTINGUE_CI"
    v2 = verdict(lo / F3, hi * F3)
    return v if v2 == v else "NAO_DISTINGUE_VIES"


def operating_characteristics(raw: dict, gates: dict) -> dict:
    """Frequencias de cada categoria de veredito por config, para as 3 regras. Limiares dos gates
    calibrados nas proprias configs nulas (0,1,12) => taxas nessas configs sao IN-SAMPLE (otimistas)."""
    thr = gates.get("thresholds_95th_pct_of_null")
    if not thr:
        return {}
    # F3: maior |vies mediano| do a0 conjunto nas configs modelo-consistentes (0,1,12), em fator
    biases = []
    for c in ("0", "1", "12"):
        if c in raw:
            b = [_bias(r["cond_joint"]["a0"], CONFIGS[int(c)]["a0"]) for r in raw[c]]
            b = [x for x in b if x is not None]
            if b:
                biases.append(abs(float(np.median(b))))
    F3 = float(10 ** max(biases)) if biases else 1.6
    out = dict(F3_bilateral_bias_factor=F3, stage2_F_range=[STAGE2_BIAS_F_LO, STAGE2_BIAS_F_HI], thresholds=thr, by_config={})
    for c, rs in raw.items():
        cnt = {"stage2_scalar": {}, "stage3_gated_joint": {}, "cut_low_only": {}}
        for r in rs:
            for key, cat in (("stage2_scalar", decide_stage2(r)), ("stage3_gated_joint", decide_stage3(r, thr, F3)),
                             ("cut_low_only", decide_cut_low(r, F3))):
                cnt[key][cat] = cnt[key].get(cat, 0) + 1
        out["by_config"][c] = dict(name=CONFIGS[int(c)]["name"], a0_true=CONFIGS[int(c)]["a0"], n=len(rs), counts=cnt)
    return out


def cross_config(raw: dict) -> dict:
    """Limiares EMPIRICOS dos gates calibrados no nulo modelo-consistente (configs 0 e 1, N=30203;
    so' 2x12 replicatas => limiar grosseiro) e aplicados as demais configs."""
    out = {}
    null_ids = [c for c in ("0", "1", "12") if c in raw]
    if not null_ids:
        return out
    def vals(c, key):
        return [r["cond_joint"]["recal"].get(key) for r in raw[c]]
    T0 = [v for c in null_ids for v in vals(c, "T_univ") if v is not None]
    G0 = [v for c in null_ids for v in vals(c, "gof_chi2") if v is not None]
    Z0 = [max(abs(np.array(z))) for c in null_ids for z in vals(c, "z_bins0_3") if z is not None]
    thrT, thrG, thrZ = (float(np.percentile(x, 95)) for x in (T0, G0, Z0))
    out["null_configs"] = null_ids
    out["thresholds_95th_pct_of_null"] = dict(T_univ=thrT, gof_chi2=thrG, zmax=thrZ, n_null=len(T0))
    out["rates_exceeding_empirical_threshold"] = {}
    for c, rs in raw.items():
        def rt(key, thr, f=lambda v: v):
            v = [f(r["cond_joint"]["recal"].get(key)) for r in rs if r["cond_joint"]["recal"].get(key) is not None]
            k = int(sum(x > thr for x in v))
            return dict(rate=k / len(v) if v else None, k=k, n=len(v), wilson95=wilson(k, len(v)))
        out["rates_exceeding_empirical_threshold"][c] = dict(
            T_univ=rt("T_univ", thrT), gof_chi2=rt("gof_chi2", thrG),
            zmax=rt("z_bins0_3", thrZ, lambda z: max(abs(np.array(z)))))
    return out


# ----------------------------------------------------------------------------------------
# Fusao de execucoes parciais (adicionado APOS as execucoes run1 (configs 0-9) e run2 (configs 10-13);
# nao altera nenhuma logica de simulacao -- ver synthetic_identifiability_v1_used_in_run1.py)
# ----------------------------------------------------------------------------------------

def part_merge(paths: list, log, result: dict) -> None:
    import hashlib
    raw: dict = {}
    sources = []
    for p in paths:
        pp = Path(p)
        b = pp.read_bytes()
        d = json.loads(b.decode("utf-8"))
        sources.append(dict(file=pp.name, sha256=hashlib.sha256(b).hexdigest(),
                            start_utc=d["meta"].get("start_utc"), end_utc=d["meta"].get("end_utc"),
                            configs=sorted(d.get("replicates", {}).keys(), key=int)))
        for k in ("unit", "explore", "retro_gof_stage2_published"):
            if k in d and k not in result:
                result[k] = d[k]
        for c, rs in d.get("replicates", {}).items():
            assert c not in raw, f"config {c} duplicada entre execucoes parciais"
            fixed = []
            for r in rs:
                r = dict(r)
                r["cond_strata"] = {int(k): v for k, v in r["cond_strata"].items()}
                fixed.append(r)
            raw[c] = fixed
        # copia do texto das partes unit/explore/retro (antes de '### PARTE main') e de todas as linhas de replicata
        lp = pp.with_suffix(".log")
        if lp.exists():
            txt = lp.read_text(encoding="utf-8").splitlines()
            head = []
            for line in txt:
                if line.startswith("### PARTE main"):
                    break
                head.append(line)
            if any(l.startswith("### PARTE unit") for l in head):
                log(f"\n##### [copiado de {lp.name}] partes unit/explore/retro #####")
                for line in head:
                    if not line.startswith("=") and not line.startswith("synthetic_identifiability.py --") and not line.startswith("inicio"):
                        log(line)
            log(f"\n##### [copiado de {lp.name}] linhas de replicata (com seeds) #####")
            for line in txt:
                if line.startswith("[") and "cfg" in line:
                    log(line)
    result["meta"]["merged_sources"] = sources
    result["summaries"] = {c: summarize(int(c), sorted(rs, key=lambda x: x["rep"])) for c, rs in sorted(raw.items(), key=lambda kv: int(kv[0]))}
    raw = {c: sorted(rs, key=lambda x: x["rep"]) for c, rs in sorted(raw.items(), key=lambda kv: int(kv[0]))}
    result["cross_config_gates"] = cross_config(raw)
    result["operating_characteristics"] = operating_characteristics(raw, result["cross_config_gates"])
    result["replicates"] = raw
    used = sorted({s for rs in raw.values() for r in rs for s in r["seeds"].values()})
    result["seeds_used_main"] = used
    log("\n### FUSAO: fontes")
    log(json.dumps(sources, indent=1, ensure_ascii=False))
    log("\n### RESUMO POR CONFIG (fusao)")
    for c in result["summaries"]:
        log(json.dumps(result["summaries"][c], indent=1, ensure_ascii=False))
    log("\n### GATES com limiar empirico (nulo = configs 0, 1 e 12)")
    log(json.dumps(result["cross_config_gates"], indent=1, ensure_ascii=False))
    log("\n### CARACTERISTICAS OPERACIONAIS das regras de decisao (contagens de veredito por config)")
    log(json.dumps(result["operating_characteristics"], indent=1, ensure_ascii=False))
    log(f"\nseeds usadas (replicatas): {len(used)} valores distintos, de {used[0]} a {used[-1]}")


# ----------------------------------------------------------------------------------------
# Execucao
# ----------------------------------------------------------------------------------------

def _worker(task):
    try:
        return run_replicate(task)
    except Exception as e:  # nao esconder falhas: registrar e propagar como resultado de erro
        import traceback
        return dict(config_id=task["config_id"], rep=task["rep"], error=repr(e), tb=traceback.format_exc())


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--part", default="all", choices=["unit", "retro", "explore", "main", "all", "smoke", "merge"])
    ap.add_argument("--workers", type=int, default=4)
    ap.add_argument("--reps", type=int, default=N_REPS_PER_CONFIG)
    ap.add_argument("--configs", default="0,1,2,3,4,5,6,7,8,9,10,11,12,13")
    ap.add_argument("--merge", nargs="+", default=None, help="JSONs de execucoes parciais a fundir (--part merge)")
    ap.add_argument("--outdir", default=str(HERE))
    ap.add_argument("--prefix", default="synthetic_identifiability")
    args = ap.parse_args()

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    log_path = outdir / f"{args.prefix}.log"
    logf = open(log_path, "w", encoding="utf-8")

    def log(msg=""):
        line = str(msg)
        print(line, flush=True)
        logf.write(line + "\n")
        logf.flush()

    log("=" * 100)
    log("synthetic_identifiability.py -- SPARC-FMULTI-STAGE3 (RASCUNHO, DISC-DEC-155). SOMENTE dados sinteticos.")
    log(f"inicio (UTC): {time.strftime('%Y-%m-%dT%H:%M:%SZ', time.gmtime())}")
    log(f"python={sys.version.split()[0]} numpy={np.__version__} scipy={__import__('scipy').__version__}")
    log(f"bloco de seeds reservado: {SEED_BLOCK_LO}-{SEED_BLOCK_HI}; "
        f"seeds de replicata: 20260955000 + 6*(10*config_id + rep) + papel(0..5) [config 0..13]; dev/unit: {DEV_SEED_BASE}-{SEED_BLOCK_HI}")
    log(f"parametros: N_MC={N_MC} N_BOOT={N_BOOT} XTOL={XTOL} PI_HIGH={PI_HIGH:.6f} bins={BIN_EDGES.tolist()} ancora={ANCHOR}")
    log("=" * 100)

    result = dict(meta=dict(script="synthetic_identifiability.py", draft="PREREGISTRATION_STAGE3_DRAFT.md",
                            n_mc=N_MC, n_boot=N_BOOT, xtol=XTOL, pi_high=PI_HIGH, a0_A=A0_A, a0_B=A0_B,
                            seed_block=[SEED_BLOCK_LO, SEED_BLOCK_HI], configs=CONFIGS,
                            n_reps_per_config=args.reps,
                            start_utc=time.strftime('%Y-%m-%dT%H:%M:%SZ', time.gmtime())))

    if args.part in ("unit", "all", "smoke"):
        log("\n### PARTE unit: testes unitarios do gerador e do pipeline")
        result["unit"] = part_unit(log)
    if args.part in ("explore", "all"):
        log("\n### PARTE explore: escolha do botao de artefato (1 replicata, N=12944, seeds dev 968-973)")
        result["explore"] = part_explore(log, args.workers)
    if args.part in ("retro", "all"):
        log("\n### PARTE retro: GOF aritmetico sobre numeros JA publicados do Estagio 2")
        result["retro_gof_stage2_published"] = part_retro_gof(log)

    if args.part == "merge":
        assert args.merge, "--merge requer caminhos"
        part_merge(args.merge, log, result)

    if args.part in ("main", "all", "smoke"):
        import multiprocessing as mp
        cfg_ids = [int(x) for x in args.configs.split(",")]
        if args.part == "smoke":
            cfg_ids = [1]
        tasks = []
        for c in cfg_ids:
            cfg = dict(CONFIGS[c])
            if args.part == "smoke":
                cfg["n_sys"] = 3000
            for rep in range(args.reps):
                tasks.append(dict(config_id=c, rep=rep, cfg=cfg,
                                  n_mc=60 if args.part == "smoke" else N_MC,
                                  n_boot=100 if args.part == "smoke" else N_BOOT))
        log(f"\n### PARTE main: {len(tasks)} replicatas ({len(cfg_ids)} configs x {args.reps}); workers={args.workers}")
        by_cfg = {c: [] for c in cfg_ids}
        t0 = time.time()
        ctx = mp.get_context("fork")
        with ctx.Pool(args.workers) as pool:
            for r in pool.imap_unordered(_worker, tasks):
                if "error" in r:
                    log(f"ERRO cfg={r['config_id']} rep={r['rep']}: {r['error']}\n{r['tb']}")
                    continue
                by_cfg[r["config_id"]].append(r)
                sc, cj = r["scalar"], r["cond_joint"]
                log(f"[{time.time() - t0:7.0f}s] cfg{r['config_id']} rep{r['rep']:02d} {r['name']}: "
                    f"f_s={sc['f_hat']:.3f} f_low={r['cond_strata'][0]['f_hat']:.3f} f_high={r['cond_strata'][1]['f_hat']:.3f} "
                    f"D={r['D_test']['D']:+.3f} a0_s={sc['a0']} a0_c={cj['a0']} "
                    f"v_s={sc['verdict']} v_c={cj['recal']['verdict']} ({r['elapsed_s']:.0f}s) seeds={r['seeds']}")
        summaries = {}
        raw = {}
        for c in cfg_ids:
            rs = sorted(by_cfg[c], key=lambda x: x["rep"])
            raw[str(c)] = rs
            if rs:
                summaries[str(c)] = summarize(c, rs)
        result["summaries"] = summaries
        result["cross_config_gates"] = cross_config(raw)
        result["operating_characteristics"] = operating_characteristics(raw, result["cross_config_gates"])
        result["replicates"] = raw
        used = sorted({s for c in cfg_ids for r in by_cfg[c] for s in r["seeds"].values()})
        result["seeds_used_main"] = used
        log("\n### RESUMO POR CONFIG")
        for c in cfg_ids:
            if str(c) in summaries:
                log(json.dumps(summaries[str(c)], indent=1, ensure_ascii=False))
        log("\n### GATES com limiar empirico (nulo = configs 0, 1 e 12)")
        log(json.dumps(result["cross_config_gates"], indent=1, ensure_ascii=False))
        log("\n### CARACTERISTICAS OPERACIONAIS das regras de decisao (contagens de veredito por config)")
        log(json.dumps(result["operating_characteristics"], indent=1, ensure_ascii=False))
        log(f"\nseeds usadas (replicatas): {len(used)} valores distintos, de {used[0]} a {used[-1]}")

    result["meta"]["end_utc"] = time.strftime('%Y-%m-%dT%H:%M:%SZ', time.gmtime())
    json_path = outdir / f"{args.prefix}.json"
    with open(json_path, "w", encoding="utf-8") as f:
        json.dump(result, f, indent=1, ensure_ascii=False, default=lambda o: o.tolist() if hasattr(o, "tolist") else str(o))
    log(f"\nJSON salvo em {json_path}")
    logf.close()


if __name__ == "__main__":
    main()
