# Pré-registro (RASCUNHO): Estágio 3 — calibração de `f_multi` condicional a RUWE (`SPARC-FMULTI-STAGE3`)

**Status:** DRAFT

**Data de criação:** 2026-10-25
**Autor (agente/sessão):** Tamesis Discovery Lab, agente de rascunho (Claude Sonnet 5.5, Claude Code), sob `DISC-DEC-155`
**Test ID:** `SPARC-FMULTI-STAGE3` (nova ETAPA de `DISC-COSMOLOGY-MOND-SPARC-004` — não é um `test_id` independente; mesma nomenclatura de `SPARC-FMULTI-STAGE1/2`).
Reabertura legítima nomeada em `01_PORTFOLIO/TEST_QUEUE.yaml:403-411` (`next_action`) e em `00_GOVERNANCE/CLAIM_LEDGER.yaml` (`DISC-CLAIM-008`, `known_gaps`).
**Autoridade para redigir:** `DISC-DEC-155` (`00_GOVERNANCE/DECISION_LEDGER.yaml:10074-10119`) — **somente rascunho**.
**Commit em que foi travado:** N/A — RASCUNHO, sem lock. Travar (LOCK), baixar o catálogo e abrir o holdout são decisões **separadas, reservadas ao autor humano** (`DISC-DEC-155`).

> **ATENÇÃO — RASCUNHO, NÃO TRAVADO, NÃO EXECUTÁVEL SOBRE DADO REAL.** Nenhum dado real foi baixado, aberto, importado ou computado na redação
> deste documento: nem `catalog.parquet`, nem Gaia, nem El-Badry, nem Hwang completo, nem `quality_filtered_sample.parquet`, nem
> `hwang_eccentricity_subset.parquet`, nem `discovery_holdout_split.json` (nem sequer as chaves do holdout). Toda a estrutura de dado citada vem de
> documentação já escrita (citada por `arquivo:linha`). Todo número novo vem de simulação **sintética** (`synthetic_identifiability.py` →
> `synthetic_identifiability.log/.json`, `analysis_summary.log`; seeds somente do bloco 20260955000-20260955999). Nenhum comando git foi usado.
> **Nenhuma alegação de Tamesis, de Problema do Milênio ou de falsificação de H_A/H_B é feita aqui.**

**Mapeamento para o template (`PREREGISTRATION_TEMPLATE.md`) e para o Estágio 2:**

| template | seção deste rascunho | observação |
|---|---|---|
| 1. Hipótese exata | §1 | física **inalterada**; muda só o modelo de contaminação |
| 2. Fonte de dado | §2 | reuso vs. re-download (documentado, não executado) |
| 3. Modelo nulo | §3 | três níveis + concorrente convencional nomeado |
| 4. Estatística de teste | §4 | procedimento passo a passo |
| 5. Critério de falsificação | §5 | regra a priori, zona de "não distingue" explícita |
| 6. Comparações múltiplas | §7 | contagem declarada antes |
| 7. O que NÃO está sendo testado | §10 | |
| (Estágio 2 §6, §8, §9, §12) | §6, §8, §9, §11 | gatilhos adversariais; linha vermelha do holdout; revalidação sintética; julgamentos abertos |
| (novo) | §0, §12 | resumo + **recomendação**; proveniência da simulação |

---

## 0. Por que esta etapa existe, o que a simulação sintética mostrou, e RECOMENDAÇÃO

### 0.1 Contexto (tudo do registro; nada recalculado)

- `SPARC-FMULTI-STAGE2` aplicou ao dado real de descoberta (30.203 sistemas) a auto-calibração escalar de `f_multi` (`f̂=0,103666`),
  obteve `a0_fit = 6,125×10⁻¹¹ m/s²` (IC95% `[4,1196; 8,5808]×10⁻¹¹`), veredito mecânico `BOTH_FALSIFIED`
  (`RESULTS_PRIMARY_STAGE2.md:78,124-125,138-155`). A reprodução independente foi `CONFIRMED` bit a bit e o debunker **não aceitou** o veredito:
  o subgrupo RUWE-alto (`max(RUWE1,RUWE2)>1,4`, 5784/30203 = 19,15%) permanece, mesmo no bin-âncora, com `+0,171` dex (IC95% `[0,126; 0,212]`) após a correção
  escalar, enquanto o RUWE-baixo fica em `-0,066` (`ADVERSARIAL_DEBUNKER_REPORT.md:179,317-335`). `DISC-DEC-032` fechou `CLOSED_INCONCLUSIVE`; `DISC-CLAIM-008` registra.
- O caminho de reabertura legítimo já escrito: "calibrar `f_multi` condicionalmente a RUWE ... exige seu proprio pre-registro" (`TEST_QUEUE.yaml:406-409`).
  `DISC-DEC-155` autoriza **só** este rascunho e o critério de refutação da frente: o rascunho é rejeitado se o referee hostil achar
  (a) caminho para tocar holdout/dado real sem lock posterior, (b) zona de veredito ambígua/reinterpretável pós-hoc, (c) validação sintética incapaz de distinguir condicional de escalar,
  (d) inconsistência com `PREREGISTRATION_STAGE2.md`/`TEST_QUEUE.yaml`.

### 0.2 O que este rascunho propõe

Manter **exatamente** a hipótese de física do Estágio 2 (H_A, H_B, `δ_obs-newt`, 5 bins, âncora = 4) e trocar o modelo de contaminação escalar por um modelo
**constante por partes em RUWE**, `f_multi,i = f_{k(i)}`, com bins pré-declarados e contagem mínima por bin (§1, §4), mais uma regra de decisão a priori com
**gates** que só podem rebaixar o veredito (§5).

### 0.3 O que a simulação sintética encontrou (14 configurações × 10 réplicas, N de 3.000 a 30.203; `analysis_summary.log` T1-T8)

Todas as afirmações abaixo são sobre **populações sintéticas** geradas pelo mesmo modelo de injeção que o pipeline usa ("inverse crime", Notas A1), exceto onde indicado.

| # | achado | números (fonte) |
|---|---|---|
| S1 | **O método condicional DISTINGUE população heterogênea de homogênea** quando o menor estrato tem ≥ ~230 sistemas no bin-âncora; falha com ~120 | detecção de `Δf=f̂₂-f̂₁` (IC95% do bootstrap do procedimento exclui 0): **10/10** com `n_âncora,alto` = 1140 (N=30.203), 487 (N=12.944), 229 (N=6.000); **2/10** com 118 (N=3.000). Falso positivo em população homogênea: **0/10** (N=30.203, cfg 0) e **0/10** (N=12.944, cfg 8) e 0/10 (`f=0`, cfg 12). O teste `D` do escalar: 10/10, 10/10, 10/10, 8/10 (T3) |
| S2 | **O `a0` do método ESCALAR NÃO é enviesado pela heterogeneidade em si** (RUWE independente de tudo): população heterogênea (`f=0,02/0,45`) e homogênea (`f=0,10`) dão o mesmo viés | viés `log10(a0_fit/a0_true)`: **+0,216 ± 0,055** (cfg 1) vs **+0,213 ± 0,028** (cfg 0), Welch p = 0,88; sem nenhuma contaminação (cfg 12): **+0,212 ± 0,035** (T1, T4). Argumento: se RUWE é independente de `g_N` e de todo o resto, a população agrupada é distribucionalmente idêntica a uma homogênea com `f̄=Σπ_k f_k` |
| S3 | **O condicional só muda o `a0` quando a composição de RUWE varia com `g_N`** (cenários hipotéticos H2) — aí o escalar erra muito e o condicional o corrige | fração RUWE-alta crescendo com `g_N`: escalar **−0,252 ± 0,029**, condicional **+0,229 ± 0,024**; decrescendo: escalar **+0,595 ± 0,027**, condicional **+0,195 ± 0,108**, corte-baixo +0,230 (T1; cfg 10/11). O veredito do Estágio 2 (regra travada) dá `BOTH_FALSIFIED` em **20/20** réplicas desses cenários **embora H_A seja verdadeiro** (T6) |
| S4 | **Todos os estimadores têm viés intrínseco de ~ +0,21 dex (fator 1,63)**, mesmo sem contaminação; o IC95% de largura ~0,16-0,18 dex **exclui o valor injetado em 0/10 réplicas** (cobertura 0/10) no N da descoberta | cfg 12 (`f=0`): escalar +0,212, condicional +0,206, corte +0,215; largura IC escalar 0,16 dex; cobertura 0/10 (T1, T2). Consequência: a regra mecânica do Estágio 2 daria `BOTH_FALSIFIED` mesmo com H_A verdadeiro, e a "camada interpretativa" (÷1,4-1,6) a rebaixa a `NAO_DISTINGUE` em 30/30 réplicas modelo-consistentes (T6): **no N da descoberta o método, tal como validado, só consegue falsificar ou declarar `NAO_DISTINGUE`** (`A_ONLY` com H_A verdadeiro só aparece em N ≤ 12.944: 2-4/10, por largura de IC) |
| S5 | **O custo do condicional é incerteza:** propagar a incerteza de `f̂` (recalibração) alarga o IC ~2× | largura IC95% (dex): escalar 0,18 / condicional-recal 0,33-0,40 / condicional-`f̂`-fixo 0,18-0,20 (N=30.203); condicional-recal 0,62 (N=12.944), 1,02 (N=6.000 — **maior que `log10(a0_B/a0_A)=0,80`: não separa H_A de H_B**), 2,10 (N=3.000) (T2, T7) |
| S6 | **Mis-especificação do estrato RUWE-alto não é corrigida pela condicionalização** (artefato aditivo de velocidade de 0,3 km/s só na verdade, ausente do mock): enviesa o `a0` conjunto *mais* que o escalar; o corte (estrato baixo isolado) fica no viés basal | cfg 3: escalar **+0,350 ± 0,032**, condicional **+0,409 ± 0,040**, corte **+0,201 ± 0,031** (p<10⁻⁵ vs cfg 1 para escalar e condicional; p=0,46 para o corte) (T1, T4); N=12.944 (cfg 13): +0,351 / +0,386 / +0,205 |
| S7 | **Os gates detectam esse artefato** (universalidade `T`, GOF), com poder 10/10 em N=30.203 e ~8/10 em N=12.944; fraco para artefato de 0,1 km/s (≤ 4/10) | limiares empíricos (percentil 95 de 30 réplicas nulas): `T*=2,45`, `χ²*=5,20`, `zmax*=1,11`. cfg 3: `T`, `χ²`, `zmax` excedem em 10/10, 10/10, 10/10; cfg 13: 7/10, 8/10, 8/10; cfg 4: 1/10, 2/10, 4/10 (`synthetic_identifiability.log`, "GATES") |
| S8 | **A regra de decisão proposta, implementada por completo (gates, gatilhos 1/2/9, CI com recalibração, largura, camada bilateral), nunca falsificou falsamente a hipótese verdadeira nas 130 réplicas em que há uma (configs 0-1, 3-13)**; a regra do Estágio 2 o fez em 10/10 (cfg 3, 10, 11), 7/10 (cfg 13), 3/10 (cfg 4)... **ao custo de poder total** | T8 (e T6). No cenário Newtoniano a regra do Estágio 2 dá `BOTH_FALSIFIED` (correto), a proposta dá `NAO_DISTINGUE[GATILHO1]` (7/10) porque o `a0` ajustado cai >1 ordem abaixo de ambos. **Com H_A verdadeiro a regra proposta completa NÃO produziu `A_ONLY` em nenhuma das 120 réplicas** (12 configs): o alargamento bilateral `F=1,65` que absorve o viés intrínseco de +0,21 dex soma +0,43 dex à largura do IC (0,3-0,4 dex com recalibração) e ultrapassa `log10(a0_B/a0_A)=0,80` (ramo L4). Só H_B verdadeiro (cfg 9) chega ao ramo final: `B_ONLY` 5/10. Sem os ramos L2/L4 (T6) haveria `A_ONLY` 4/10 (N=30.203) e 7/10 (N=12.944) — frequência governada pela largura do IC alargado, não por resolução |
| S9 | **Retrospectiva aritmética sobre números JÁ publicados** (não é teste; Notas V5): os gates reprovariam no dado de descoberta | universalidade a partir da tabela do debunker (`ADVERSARIAL_DEBUNKER_REPORT.md:319-323`): `Δ_b = [0,229; 0,219; 0,187; 0,102]`, `z=[4,6; 4,7; 4,2; 2,3]`, `T=66,2`; GOF da forma AQUAL sobre o `δ` corrigido do Estágio 2 (`RESULTS_PRIMARY_STAGE2.md:95-99`): `χ²=45,6` no `a0_fit` publicado (`synthetic_identifiability.log`, `RETRO-*`). O jitter de 0,3 km/s sintético dá `Δ_b=[0,29; 0,19; 0,09; 0,04]` — **qualitativamente** semelhante, sem provar causa |

**Resposta à hipótese da frente (`DISC-DEC-155`):** "é possível pré-declarar um `f_multi(RUWE)` cuja recuperação sintética distinga o caso condicional do escalar?" —
**Em parte.** (i) *Distinção estatística de heterogeneidade*: sim, com poder ≥ 0,8 e falso positivo ≤ 0,1 para `n_âncora,alto ≥ ~230` (S1). (ii) *Distinção em termos de `a0` recuperado*: **só** quando a fração de RUWE-alto
varia com `g_N` (S3) — fato do dado real que **não conheço** — e **não** quando a heterogeneidade é pura (S2) nem sob mis-especificação (S6). (iii) O viés intrínseco (S4) limita qualquer das duas versões a falsificar ou declarar `NAO_DISTINGUE`.
O critério (c) de refutação da frente ("validação sintética incapaz de distinguir condicional de escalar") **aplica-se parcialmente**: aplica-se ao caso-base (S2), não ao caso H2 (S3).

### 0.4 RECOMENDAÇÃO (não vinculante; decisão do autor humano — Seção 11)

**Não travar este rascunho como está** (`Status` permanece `DRAFT`). Razões, em ordem de peso:

1. **Viés intrínseco do estimador de `a0` (S4, J8):** +0,21 dex independente de contaminação; enquanto existir, nenhum veredito de *suporte* a H_A ou H_B é alcançável e a regra do Estágio 2 produz falsificação espúria de H_A verdadeiro quando há dependência com `g_N` (S3, S8).
   Isto pede uma revisão do **estimador** (J8), que está fora do escopo "calibração condicional de `f_multi`" — logo, é matéria de decisão humana e de novo desenho, não de ajuste neste documento.
2. **A premissa motivadora não se traduz em viés de `a0` no caso-base (S2).** A evidência do debunker (excesso no âncora do RUWE-alto) é *exatamente* o que `f` heterogêneo independente produz (cfg 1 reproduz `D≈+0,28`), sem afetar o `a0` agrupado. Só a dependência de `g_N` (S3) ou a forma da contaminação (S6) o afetariam — e o que o debunker mostra de *forma* (S9) aponta para S6, que o modelo condicional **não** resolve.
3. **Execução na amostra de descoberta provavelmente resultaria em `NAO_DISTINGUE`** (S9), e a forma/borda do estrato foram motivadas por essa mesma amostra (J1: dupla utilização). O valor informativo esperado de uma execução na descoberta é baixo.
4. **Poder:** com o N do holdout (12.944) o IC conjunto (0,62 dex) mal cabe em 0,80 dex; com ~6.000 não cabe (S5). Isso informa que um confirmatório só no holdout tem margem estreita. E, com a regra completa (T8), H_A verdadeiro **não** alcança `A_ONLY` em nenhum N testado (S8).
5. **19 julgamentos abertos (Seção 11)**, vários estruturais (J1, J2, J8), e **10 inconsistências com o registro** (+1 incidente de procedimento meu; `STAGE3_DRAFT_NOTES.md`, §5), incluindo a lacuna de "holdout lido em memória" (I-5).

**Se o autor decidir prosseguir**, a ordem sugerida é: J8 (estimador) → J1 (papel da descoberta) → J2 (estimador primário; o **corte** é o robusto em todos os cenários sintéticos) → J4-J6 → só então revalidar (Seção 9) e considerar o LOCK.
Uma versão **só-diagnóstica** (J17b: caracterizar o estrato RUWE-alto com os gates de universalidade/GOF, sem veredito H_A/H_B) é a que a simulação melhor sustenta, mas também é decisão do autor.

---

## 1. Hipótese exata

### 1.1 A hipótese de física NÃO é reformulada (idêntica a SPARC-002/003/004 e Estágio 2)

- **H_A:** o `a0` que melhor explica `δ_obs-newt(bin)` corrigido por multiplicidade (§4) é compatível, dentro do IC95%, com
  `a0^A = cH0/(2π) ≈ 1,082288×10⁻¹⁰ m/s²` ("Ponte Holográfica") — `PREREGISTRATION_STAGE2.md:64-67`, `PREREGISTRATION.md:65-71`.
- **H_B:** o mesmo `a0` é compatível com `a0^B = cH0 ≈ 6,800218×10⁻¹⁰ m/s²` ("MOND Emergence") — `PREREGISTRATION_STAGE2.md:68-70`, `PREREGISTRATION.md:72`.
- Mesmos valores exatos, mesma função `ν(x)=1/(1-e^{-√x})`, mesmos cortes de qualidade, mesmo split discovery/holdout, mesmas bordas de bin e bin-âncora (§2).
- **Nota de identificabilidade honesta (carregada sem diluir):** `a0^A=cH0/2π` reproduz uma coincidência numérica já conhecida na literatura MOND (Milgrom, arXiv:2001.09729;
  `PREREGISTRATION.md:78-84`). H_A sobreviver não é evidência específica de Tamesis.
- **Discriminating observable (`METHODOLOGY_EXTENSIONS.md` §1):** a **universalidade de `a0` entre estratos de RUWE** depois do âncora por estrato (§4.3 E5). Modelo concorrente nomeado que ele separa de gravidade:
  `M_conv` = companheiras ocultas **mais** artefato astrométrico dependente de RUWE não modelado (§3).

### 1.2 O que muda: o modelo de contaminação `f_multi(RUWE)`, em forma paramétrica explícita

```text
R_i       = max(RUWE1_i, RUWE2_i)
k(i)      = 1  se R_i <= 1,4        (estrato RUWE-baixo)
          = 2  se R_i  > 1,4        (estrato RUWE-alto)                 [K = 2, primário; J3]
f_multi,i = f_{k(i)},   (f_1, f_2) in [0 ; 0,9]^2                        [função constante por partes em R; 2 parâmetros livres]
```

O restante do modelo de injeção é **idêntico** ao do Estágio 2: Eqs. 11-13 de Chae (`PROVENANCE_CHAE_EQS.md:47-90`), componentes afetados 40/30/30%, `γ_M=-0,7`, wobble de fotocentro simplificado, desprojeção, ruído astrométrico simétrico (`PREREGISTRATION.md` §5b).
**Contagem mínima por bin:** `n_{k,4} ≥ N_A_MIN` e `n_k ≥ N_S_MIN` (valores propostos na §5.1; regra de parada se violada). **Por que esta forma e não outra:** §4.2.

### 1.3 O que a frente afirma, decomposto em enunciados testáveis sinteticamente (DISC-DEC-155)

- **E-S1:** o procedimento condicional distingue uma população com `f` dependente de RUWE de uma homogênea (poder e falso positivo). → **sim**, para `n_âncora,alto ≳ 230` (§0.3 S1).
- **E-S2:** o `a0` escalar é enviesado na população heterogênea e o condicional recupera o injetado. → **não, no caso-base** (S2); **sim**, quando a composição varia com `g_N` (S3).
- **E-S3:** o procedimento é identificável sob contaminação não modelada. → **não** (S6); os gates o **detectam** (S7) mas não o **corrigem**.

Nenhum destes enunciados é uma afirmação sobre o mundo físico — são afirmações sobre o comportamento do **procedimento** em populações geradas por construção.

## 2. Fonte de dado (nada é baixado, aberto ou computado por este rascunho)

Esta seção é **documentação**, não execução. Todas as colunas e contagens vêm de
arquivos de texto já escritos (citados), nunca de uma leitura do dado.

### 2.1 O que seria REAUTILIZADO (já baixado e commitado por Estágios anteriores; nenhum novo download)

| item | caminho (relativo a `02_TESTS/`) | uso no Estágio 3 | fonte da informação |
|---|---|---|---|
| `quality_filtered_sample.parquet` | `COSMOLOGY_WIDE_BINARIES/data/` | `M1_Msun`, `M2_Msun`, erros de PM por componente (já usados) e, **novo uso**, `RUWE1`, `RUWE2` como variável de estratificação (no Estágio 2 eram só "checagem de consistência") | `TEST_QUEUE.yaml:20`; `fmulti_selfcal_stage1/METHODOLOGY_ADDENDUM.md:95-100` (tabela de insumos) |
| `hwang_eccentricity_subset.parquet` | `COSMOLOGY_WIDE_BINARIES_MC_DEPROJECTION/data/` | `e`, `e0`, `e1`, `alpha`, `dpm_sig` (Gap (a), inalterado) | `PREREGISTRATION.md:105-119`; `METHODOLOGY_ADDENDUM.md:99` |
| `discovery_holdout_split.json` | `COSMOLOGY_WIDE_BINARIES/data/` | **somente** as chaves `discovery_pair_ids` e `n_discovery` (ver Seção 8) | `PREREGISTRATION.md:98-102`; `PREREGISTRATION_STAGE2.md:441-470` |
| `mamajek_mass_luminosity.tsv` | `COSMOLOGY_WIDE_BINARIES/data/` | relação massa-magnitude (inalterada) | `PREREGISTRATION.md:103-104` |
| módulos travados `selfcal_pipeline.py`, `companion_injection.py`, `build_synthetic_population.py` | `fmulti_selfcal_stage1/analysis/` | pipeline de auto-calibração, **aplicado por estrato** (sem edição: `calibrate_f_multi` já funciona sobre qualquer subconjunto de sistemas) | `PREREGISTRATION_STAGE2.md:179-222` |
| `deprojection_common.py`, `delta_obs_newt.py` (só `BIN_EDGES_LOG_GN_SPARC003`) | `COSMOLOGY_WIDE_BINARIES_MC_DEPROJECTION/analysis/` | desprojeção e bordas de bin (inalteradas) | `PREREGISTRATION_STAGE2.md:118-133` |

Amostra: os **30.203** sistemas de descoberta (`PREREGISTRATION_STAGE2.md:87`; `RESULTS_PRIMARY_STAGE2.md:25`),
mesmos cortes de qualidade, mesmas bordas de bin `[-11,7012; -9,1728; -8,4667; -7,9752; -7,5548; -6,5224]`
(`PREREGISTRATION_STAGE2.md:124`) e mesmo bin-âncora = 4 (`PREREGISTRATION_STAGE2.md:134-138`).
**Nenhum registro novo de dado.** O único insumo conceitualmente novo é o uso de `RUWE1`/`RUWE2` como estratificador.

### 2.2 O que TERIA de ser baixado/aberto de novo se o autor humano escolher certas opções (documentar, NÃO baixar)

| necessidade | quando surge | por que não está disponível hoje | custo/risco |
|---|---|---|---|
| Colunas Gaia de qualidade astrométrica adicionais (`astrometric_excess_noise`, `ipd_*`, `visibility_periods_used`, ...) como estratificador alternativo ou contínuo | só se o autor rejeitar RUWE como único estratificador (J3/J15, Seção 11) | **não verificado se existem** em `quality_filtered_sample.parquet`: as únicas colunas documentadas são as de `METHODOLOGY_ADDENDUM.md:95-100`; não abri o arquivo | exigiria nova consulta Gaia/El-Badry + nova proveniência (`AGENTS.md` passo 5) |
| `catalog.parquet` bruto (El-Badry, ~1,8 M de pares pré-corte de qualidade) para obter a distribuição de RUWE **antes** dos cortes | só se o autor quiser o baseline populacional de RUWE (J16) | o debunker deliberadamente **não** o abriu (`ADVERSARIAL_DEBUNKER_REPORT.md:94-100,470-474`); contém fisicamente os mesmos sistemas do holdout, sem rótulo | exige autorização humana separada (mesma lógica da Seção 8) |
| Catálogo externo de multiplicidade de alta ordem / velocidade radial para validar `f_k` contra um gabarito independente | só para validação externa (opcional) | "NÃO disponível/não buscado" (`METHODOLOGY_ADDENDUM.md:101`) | fora do escopo |
| Extrato "só-descoberta" da tabela de qualidade (arquivo novo com apenas as 30.203 linhas) | se o autor adotar J14 | **não existe**; criá-lo exige ler a tabela completa + a lista de descoberta uma vez | passo de preparação de dado com autorização própria; produz sha256 |
| Linhas do holdout (12.944) | **nunca neste Estágio** | selado | Seção 8 |

**Fato do registro que este rascunho precisa tornar explícito (ver `STAGE3_DRAFT_NOTES.md`, inconsistência I-5):**
`run_stage2_primary_analysis.py:178` carrega `quality_filtered_sample.parquet` inteiro com `pd.read_parquet`
e só depois filtra por `discovery_pair_ids` (linha 200). Pela documentação (`PREREGISTRATION.md:94-97`:
"43.147 sistemas pós-corte"; `PREREGISTRATION.md:98-100`: 30.203 + 12.944), o arquivo contém também as linhas do holdout.
Logo, a frase "nenhum dado físico de nenhum sistema do holdout foi lido" (`RESULTS_PRIMARY_STAGE2.md:246-249`)
vale no sentido "nenhuma coluna do holdout foi *computada ou usada*", não no sentido "nenhuma linha foi carregada em memória".
Não verifiquei o conteúdo do parquet (proibido); isto é inferência a partir da documentação.

## 3. Modelo nulo / hipótese concorrente

Três níveis, para não confundir o que é testado:

1. **Nulo de física (inalterado, `PREREGISTRATION_STAGE2.md:140-152`):** consistência interna entre H_A e H_B contra um canal
   físico independente; sob gravidade Newtoniana + contaminação *corretamente modelada*, `δ_obs-newt(bin) = 0` em todos os bins
   depois da correção por multiplicidade. Não é um teste "MOND vs. Newton" no sentido do Artigo A de Chae.
2. **Nulo da frente (o que o Estágio 3 muda): `f_multi` homogêneo na população** — o modelo escalar do Estágio 2.
   O concorrente é o modelo condicional `f_multi,i = f_{k(i)}` da Seção 1. A pergunta da frente é se os dois são
   **distinguíveis** com o poder disponível, e se a distinção muda o `a0` recuperado (Seção 0, resultados sintéticos).
3. **Concorrente convencional nomeado (`METHODOLOGY_EXTENSIONS.md` §1 e §5):** `M_conv` = companheiras ocultas **mais** um artefato
   astrométrico dependente de RUWE que o modelo de injeção (massa + wobble de fotocentro, Eqs. 11-13 de Chae) **não** representa
   (solução astrométrica degradada, erro de PM subestimado em RUWE alto — `ADVERSARIAL_DEBUNKER_REPORT.md:164-173,350-358`).
   `M_conv` prediz uma assinatura específica e testável: o excesso do estrato RUWE-alto **não** se anula nos bins de menor `g_N`
   depois de ancorado, i.e. viola a universalidade de `a0` entre estratos (§4.3, passo E5). Gravidade (qualquer `a0`) prediz o oposto.
   A simulação (Seção 0, configs 3/4/13) mostra que um artefato aditivo de velocidade, não modelado, reproduz qualitativamente essa assinatura.

## 4. Estatística de teste e procedimento exato (PROPOSTO; parâmetros marcados «J» são julgamentos abertos da Seção 11)

Nada abaixo foi executado sobre dado real. O protótipo sintético é `synthetic_identifiability.py` (funções `run_replicate`,
`fpb`, `fit_a0_joint`); a implementação real deve ser reescrita por agente separado e comparada ao protótipo (Seção 9).

### 4.1 Definições

- `R_i = max(RUWE1_i, RUWE2_i)` (mesma definição do debunker, `ADVERSARIAL_DEBUNKER_REPORT.md:179` e `RESULTS_PRIMARY_STAGE2.md:225-228`).
- **Estratos (K=2, primário, J3):** `k(i)=1` se `R_i ≤ 1,4`; `k(i)=2` se `R_i > 1,4`. A borda 1,4 é o limiar convencional de "ajuste
  astrométrico ruim" já usado como `RUWE_FLAG_THRESHOLD = 1.4` no Estágio 2 (`run_stage2_primary_analysis.py:105`) e na checagem do debunker.
- **Modelo de contaminação (a forma paramétrica explícita):** `f_multi,i = f_{k(i)}`, função **constante por partes** em `R`, com
  `K=2` parâmetros livres `(f_1, f_2) ∈ [0, 0,9]²`. O restante do modelo de injeção (Eqs. 11-13 de Chae, 40/30/30% de componentes
  afetados, `γ_M = -0,7`, wobble simplificado) é **idêntico** ao do Estágio 2.
- Bins de `g_N`, bin-âncora (=4), `δ_obs-newt`, `ν` "simple", cortes de qualidade: **idênticos** ao Estágio 2 (Seção 2).

### 4.2 Por que esta forma (e não outra)

1. **Menos parâmetros que a evidência suporta.** A calibração de Chae usa **um** número por população (o bin-âncora). Cada estrato
   adicional custa um bin-âncora com menos sistemas; a simulação (Seção 0, T7) mede o preço em N. Constante por partes com K=2 é a
   forma mais simples que permite a distinção que a evidência existente pede (subgrupo RUWE-alto vs. resto).
2. **Mesmo mecanismo de calibração.** Cada `f_k` é calibrado pela *mesma* regra verbatim de Chae
   (`PROVENANCE_CHAE_EQS.md:115-167`: iterar até `δ_obs-newt = 0` no bin de maior aceleração), apenas **dentro do estrato**.
   Isto reaproveita `calibrate_f_multi()` sem edição (Seção 2).
3. **RUWE é um proxy ruidoso, não a causa.** A relação entre RUWE e fração de companheiras ocultas é não linear e desconhecida;
   uma forma contínua (p.ex. logística em `log(RUWE-1)`) teria parâmetros que o único bin-âncora **não** consegue identificar.
4. **A forma não é nova para o registro:** é exatamente a divisão que o debunker aplicou (`ADVERSARIAL_DEBUNKER_REPORT.md:304-335`).
   **Isto é também a sua fraqueza:** a forma foi sugerida olhando a amostra de descoberta (J1, Seção 11) e, portanto, não pode
   ser tratada como hipótese "cega" nessa amostra.
5. **Resultado sintético que limita a justificativa (Seção 0):** se `f` depende de RUWE mas RUWE é independente de `g_N` e de todo o
   resto, o estimador **escalar** de `a0` não é enviesado (a população agrupada é distribucionalmente idêntica a uma homogênea com
   `f̄ = Σ π_k f_k`). A forma condicional só muda o `a0` quando RUWE carrega informação sobre `g_N` (ou sobre a *forma* da contaminação).
   Logo, a motivação de "o escalar está mal especificado" **não** se traduz automaticamente em "o `a0` escalar está enviesado".

**Contagem mínima por bin (pré-declarada, a preencher na Seção 5.1 com o resultado da simulação):** `n_{k,4} ≥ N_A_MIN` (sistemas do
estrato `k` no bin-âncora) e `n_k ≥ N_S_MIN`; se violada para algum estrato, o Estágio 3 **não é executado** (`NAO_EXECUTAVEL`),
nunca se funde estratos *depois* de ver as contagens (regra de fusão, se houver, é pré-declarada e determinística — J4).

### 4.3 Passos

```text
E0  PORTÃO DE CONTAGEM (primeiro e único acesso a dado real antes do LOCK dos resultados): contagens n_k e n_{k,b} (k=1..K, b=0..4)
    da AMOSTRA DE DESCOBERTA já filtrada pelos 30.203 ids. Nada além de contagens. Se n_{k,4} < N_A_MIN ou n_k < N_S_MIN: PARAR.
E1  CALIBRAÇÃO POR ESTRATO:  calib_k = calibrate_f_multi(subconjunto k, anchor_bin=4, f_lo=0, f_hi=0,9, xtol, n_mc, seed_k,
    include_wobble=True, return_raw=True)  ->  f̂_k, converged_k, δ_{k,b}, gN_{k,b} (b=0..4).
E2  HETEROGENEIDADE (informativo):  Δf = f̂_2 - f̂_1  e  D = δ_{2,âncora}(f̂_escalar) - δ_{1,âncora}(f̂_escalar), com ICs do bootstrap (E4).
E3  BOOTSTRAP DO PROCEDIMENTO (FPB): para cada k, uma segunda execução pareada em f̂_k ± 0,15 (mesma seed_k+777 => sorteios aninhados);
    em cada réplica b=1..B: reamostrar sistemas e realizações MC; RECALIBRAR f_k por interpolação linear em f entre as duas
    execuções (clip em [0, 0,9]); recomputar δ_{k,b}. (Erro da interpolação linear medido em U3: <= 0,004 dex.)
E4  a0 CONJUNTO:  â0 = argmin Σ_{k,b} w_{k,b} [δ_{k,b} - log10 ν(gN_{k,b}/a0)]²,  w_{k,b} = n_{k,b} (pesos determinísticos, J10);
    IC95% = percentis 2,5/97,5 dos a0 refeitos em cada réplica do FPB.
E5  GATE DE UNIVERSALIDADE (a0-livre):  Δ_b = δ_{2,b} - δ_{1,b}, b=0..3;  z_b = Δ_b / sd_FPB(Δ_b);  T = Σ_b z_b²;  zmax = max_b |z_b|.
E6  GATE DE QUALIDADE DO AJUSTE:  χ² = Σ_{k,b<=3} [(δ_{k,b} - modelo)/sd_FPB(δ_{k,b})]².
E7  ALTERNATIVAS DESCRITIVAS (não entram na decisão, salvo via gatilho 9, Seção 6):  â0 só do estrato baixo (corte), â0 por estrato.
E8  Regra de decisão da Seção 5.
```

O ajuste de `a0` por estrato e o do estrato baixo isolado usam `fit_a0()` e `bootstrap_a0_refit()` do Estágio 2 sem alteração.

### 4.4 Parâmetros numéricos

| parâmetro | valor da simulação | valor proposto para o dado real | status |
|---|---|---|---|
| `n_mc` | 100 | 200 (como `PREREGISTRATION_STAGE2.md:306`) | J9 |
| réplicas de bootstrap `B` | 300 | 2000 (como `PREREGISTRATION_STAGE2.md:307`) | J9 |
| `xtol` da bisseção | `2×10⁻³` | `5×10⁻⁴` (default, `PREREGISTRATION_STAGE2.md:303`) | J9 |
| passo de f do FPB | 0,15 | 0,15 | — |
| `include_wobble` | True | True (inalterado) | — |
| `real_gets_astrometric_noise` | False (n/a na sim.) | False — **não sobrescrever** (`PREREGISTRATION_STAGE2.md:207-217`) | — |
| seeds do Estágio 3 real | n/a | **a definir pelo autor**, fora do bloco 20260955000-999 (reservado a simulação) | J9 |
| `N_A_MIN`, `N_S_MIN` | — | Seção 5.1 | J4 |
| limiares `T*`, `zmax*`, `χ²*` | 95º percentil dos nulos sintéticos (grosseiro) | recalibrar com ≥ 200 réplicas nulas no checklist (Seção 9) | J6 |
| fator de viés `F` | Seção 5 | Seção 5 | J7/J8 |

## 5. Critério de decisão a priori (falsificação / suporte / zona de "não distingue")

Princípio herdado de `PREREGISTRATION_STAGE2.md:333-353` e `PREREGISTRATION.md` §7e: **toda camada adicional só pode REBAIXAR o veredito para `NAO_DISTINGUE`;
nenhuma abre caminho novo para aceitar H_A ou H_B.** Não é permitido reinterpretar um resultado ambíguo como suporte depois do fato.

### 5.1 Pré-condições de executabilidade e N mínimo por bin (sintético; valores PROPOSTOS, J4)

Critérios de identificabilidade declarados no script **antes** da execução final (`summarize`, bloco `criteria`):
**C1** — detectar heterogeneidade com poder ≥ 0,80 (`D` e `Δf`) e, em população homogênea, falso positivo ≤ 0,10; **C2** — |viés mediano do `a0` conjunto| ≤ 0,30 dex (a mesma tolerância de
`validate_b_recover_a0_with_contamination.py:41`); **C3** — largura mediana do IC95% conjunto < `log10(a0_B/a0_A)` = **0,80 dex** (senão o IC nem em princípio separa H_A de H_B).

| N total | `n_âncora` do estrato alto (média; mín.) | detecção `Δf` / `D` | viés `a0` escalar / condicional (dex) | largura IC cond.-recal (dex, média) | C1 | C2 | C3 |
|---|---|---|---|---|---|---|---|
| 3.000 (cfg 7) | 118 (106) | **2/10** / 8/10 | +0,31 ± 0,16 / +0,29 ± 0,14 | 2,10 | **falha** | **falha** (mediana +0,34) | **falha** |
| 6.000 (cfg 6) | 229 (206) | 10/10 / 10/10 | +0,27 ± 0,09 / +0,25 ± 0,08 | 1,02 | passa | passa | **falha** |
| 12.944 (cfg 5) | 487 (457) | 10/10 / 10/10 | +0,22 ± 0,07 / +0,21 ± 0,08 | 0,62 | passa | passa | passa |
| 30.203 (cfg 1) | 1.140 (1.115) | 10/10 / 10/10 | +0,22 ± 0,06 / +0,20 ± 0,06 | 0,40 | passa | passa | passa |

(Fonte: `analysis_summary.log` T7 e `synthetic_identifiability.log`, `criteria`. Os "N" são os do catálogo sintético com contagens iguais por bin de `g_N`; a fração do estrato alto é 0,1915.)

- **Menor configuração testada que satisfaz C1-C3:** `n_âncora,alto ≈ 487` (mínimo observado 457), N = 12.944 (cfg 5). **Abaixo de ~230 reprova C3; abaixo de ~120 reprova também C1 e C2.**
  Não testei valores intermediários entre 229 e 457: o limiar real está nesse intervalo.
- **Proposta provisória (J4):** `N_A_MIN = 450` sistemas no bin-âncora do menor estrato **e** `n_{k,b} ≥ 450` em toda célula `(k,b)` usada no ajuste; `N_S_MIN = 5×450 = 2.250` no menor estrato.
  Se violado ⇒ `NAO_EXECUTAVEL` (nunca fundir estratos depois de ver as contagens; regra de fusão, se houver, é pré-declarada — J4).
- **Em termos da amostra real documentada:** com `π_alto=0,1915`, a descoberta (30.203) esperaria ~1.150 por célula do estrato alto — folga de ~2,5×; um holdout de 12.944 esperaria ~490 — **sem folga** (mín. 457 em 10 sorteios sintéticos).
  Se a fração RUWE-alta variar com `g_N`, as células de canto encolhem: cfg 10 (fração crescente) tem o estrato alto com `[240; 491; 968; 1.628; 2.483]` sistemas nos bins 0-4 (mín. 203 no bin 0);
  cfg 11 (decrescente) tem `[2.494; 1.628; 970; 469; 236]` (mín. 213 no **bin-âncora**), e o IC conjunto cresce para 2,56 dex de largura (cfg 11; T2) — **a contagem por célula, não a total, decide.**
- **Pré-condições adicionais:** checklist da §9 integralmente passado; entrada de dado = extrato só-descoberta (ou garantia mais fraca declarada, J14); `fit_a0` com guarda `a0>0`.

### 5.2 Regra de decisão em camadas (função determinística; cada ramo enumerado)

Seja `[lo, hi]` o IC95% (percentis 2,5/97,5 do bootstrap do procedimento com recalibração, §4.3 E3-E4) de `â0` conjunto. Defina
`S_A = (lo ≤ a0_A ≤ hi)`, `S_B = (lo ≤ a0_B ≤ hi)` e a classe mecânica `M`:
`A_ONLY = S_A ∧ ¬S_B`; `B_ONLY = ¬S_A ∧ S_B`; `BOTH_SURVIVE = S_A ∧ S_B`; `BOTH_FALSIFIED = ¬S_A ∧ ¬S_B`.
Seja `F` o fator de viés bilateral e `M'` a mesma classe calculada sobre `[lo/F, hi·F]`. Seja `W' = log10(hi·F / (lo/F))`.

```text
V(dados) =
  L0  NAO_EXECUTAVEL      se qualquer pré-condição da §5.1 / §9 falha                                  (nenhum cálculo sobre dado real prossegue)
  L1  NAO_DISTINGUE[GATE] se  T > T*  ou  zmax > zmax*  ou  chi2 > chi2*                               (§4.3 E5, E6)
  L2  NAO_DISTINGUE[GATILHO k]  se qualquer gatilho 1, 2(estrato alto, ou estrato baixo na borda superior), 6, 7, 9 ou 10 da §6 dispara
  L3  NAO_DISTINGUE[CI]   se  M = BOTH_SURVIVE
  L4  NAO_DISTINGUE[LARGURA] se  W' >= log10(a0_B/a0_A) = 0,798 dex                                     (o IC alargado não separa H_A de H_B nem em princípio)
  L5  NAO_DISTINGUE[VIES] se  M' != M                                                                   (a camada de viés resgata A e/ou B, ou as mistura)
  L6  caso contrário  V = M  ∈ {A_ONLY, B_ONLY, BOTH_FALSIFIED}                                         (SEMPRE sujeito à reexecução adversarial + debunker, §6)
```

`F` proposto: **`F = 1,65`** (= `10^max|viés mediano do a0 conjunto|` nas configs sintéticas modelo-consistentes 0, 1 e 12; valor calculado em `operating_characteristics`: `F3 = 1,653`), **bilateral** — diferente do Estágio 2, que só dividia (÷1,4-1,6)
e portanto não protegia contra viés negativo (cfg 10: −0,25 dex ⇒ `BOTH_FALSIFIED` em 10/10 réplicas). `T* = 2,45`, `zmax* = 1,11`, `chi2* = 5,20` são os percentis 95 de 30 réplicas nulas sintéticas (grosseiros; **recalibrar com ≥ 200 réplicas**, §9.2, J6).
Estes valores são **provisórios**; o que é pré-declarado é a *estrutura* e a obrigação de fixá-los antes do LOCK.

### 5.3 Zona de "não distingue" — explícita e exaustiva

`V ∈ {NAO_EXECUTAVEL, NAO_DISTINGUE[GATE], NAO_DISTINGUE[GATILHO k], NAO_DISTINGUE[CI], NAO_DISTINGUE[LARGURA], NAO_DISTINGUE[VIES]}` — isto é, **tudo** que não seja o ramo L6.
Nenhum desses resultados pode ser relatado, citado ou arquivado como suporte a H_A, a H_B, a Tamesis ou a "falsificação parcial". Também **não** é permitido:
(i) mudar K, a borda de RUWE, o estimador, os pesos, `F`, os limiares ou o conjunto de gates depois de ver qualquer resultado real; (ii) escolher entre as variantes de sensibilidade a que fornece o veredito mais favorável — qualquer mudança de categoria entre variantes rebaixa a `NAO_DISTINGUE` (§6 gatilho 10);
(iii) executar mais de uma vez sobre a descoberta (§7); (iv) relatar o estimador "corte" ou o `a0` por estrato como resultado primário.
Se os gates tiverem poder baixo no N real (checklist 9.5), isso é **documentado** como limitação, não corrigido depois.

### 5.4 Linguagem permitida por veredito

- `A_ONLY`: "H_B é excluída por este canal, sob este modelo de contaminação; H_A não é excluída (consistência, não confirmação)".
- `B_ONLY`: simétrico. `BOTH_FALSIFIED`: "ambas excluídas por este canal, sob este modelo — sujeito ao debunker". `NAO_DISTINGUE[*]`: "este canal não distingue H_A de H_B".
- **Proibido:** "Tamesis confirmado/refutado/favorecido", "H_A detectada", qualquer menção a Problema do Milênio (`AGENTS.md`, Proibições).

### 5.5 Características operacionais da regra (sintético; `analysis_summary.log` T8 e T6; `synthetic_identifiability.log` "CARACTERISTICAS OPERACIONAIS")

A tabela usa a regra **implementada o mais completamente possível no pós-processamento** (T8): gates L1, gatilhos 1, 2 (estrato alto na borda, ou estrato baixo na borda superior) e 9, ramos L3-L5 com `F=1,65`.
**Não implementados nas simulações:** L0 (pré-condições), gatilho 6 (sensibilidade a `B`/seed), gatilho 10 (variantes de sensibilidade), gatilhos 3 e 8 (só exigem escrutínio) e o gatilho 7 (dispararia por construção nas configs homogêneas 0, 8, 12).
Esses ramos só podem **rebaixar**; logo, a coluna "proposto" é um **limite superior** das conclusões alcançáveis. Limiares dos gates calibrados nas próprias configs nulas (0, 1, 12) ⇒ taxas nelas são **in-sample**. 10 réplicas por config.
Entre colchetes, o ramo de `NAO_DISTINGUE`; "(falsa)" = falsificação da hipótese que era verdadeira.

| config (verdade injetada) | Estágio 2 (escalar + camada ÷1,4-1,6) | **Estágio 3 proposto** (L1-L5, `F=1,65` bilateral) |
|---|---|---|
| 0 homogêneo, H_A (N=30.203) | `[VIES]` 10 | `[VIES]` 9, `[GATE]` 1 |
| 1 heterogêneo, H_A | `[VIES]` 10 | `[LARGURA]` 8, `[GATE]` 1, `[VIES]` 1 |
| 12 sem contaminação, H_A | `[VIES]` 10 | `[VIES]` 6, `[GATILHO2]` 3, `[GATE]` 1 |
| 2 heterogêneo, **Newtoniano** | `BOTH_FALSIFIED` 10 | `[GATILHO1]` 7, `[GATE]` 2, `[LARGURA]` 1 |
| 9 heterogêneo, **H_B** | `[VIES]` 10 | **`B_ONLY` 5**, `[VIES]` 5 |
| 10 composição cresce com `g_N`, H_A | **`BOTH_FALSIFIED` 10 (falsa)** | `[VIES]` 8, `[GATE]` 2 |
| 11 composição decresce com `g_N`, H_A | **`BOTH_FALSIFIED` 10 (falsa)** | `[CI]` 4, `[GATE]` 4, `[LARGURA]` 2 |
| 3 artefato 0,3 km/s, H_A | **`BOTH_FALSIFIED` 10 (falsa)** | `[GATE]` 10 |
| 4 artefato 0,1 km/s, H_A | `[VIES]` 7, **`BOTH_FALSIFIED` 3 (falsa)** | `[GATE]` 5, `[LARGURA]` 4, `[VIES]` 1 |
| 13 artefato 0,3 km/s, H_A (N=12.944) | **`BOTH_FALSIFIED` 7 (falsa)**, `[VIES]` 3 | `[GATE]` 9, `[LARGURA]` 1 |
| 5 heterogêneo, H_A (N=12.944) | `A_ONLY` 2, `[VIES]` 8 | `[LARGURA]` 9, `[GATE]` 1 |
| 8 homogêneo, H_A (N=12.944) | `A_ONLY` 2, `[VIES]` 8 | `[LARGURA]` 7, `[GATE]` 3 |
| 6 heterogêneo, H_A (N=6.000) | `[VIES]` 7, `A_ONLY` 2, **`BOTH_FALSIFIED` 1 (falsa)** | `[LARGURA]` 8, `[GATE]` 2 |
| 7 heterogêneo, H_A (N=3.000) | `A_ONLY` 4, `[VIES]` 4, **`BOTH_FALSIFIED` 2 (falsa)** | `[CI]` 6, `[LARGURA]` 3, `[GATE]` 1 |

Leitura honesta:
(a) **Segurança:** a regra proposta não produziu falsificação falsa (130 réplicas com hipótese verdadeira) nem suporte espúrio com a verdade trocada ou Newtoniana; a regra do Estágio 2 produziu falsificação falsa em cfg 3/10/11 (10/10 cada), cfg 13 (7/10), cfg 4 (3/10), cfg 7 (2/10), cfg 6 (1/10).
(b) **Poder:** com H_A verdadeiro, a regra proposta completa **nunca** chega a `A_ONLY` (0/120). Só H_B verdadeiro (a0 maior ⇒ IC mais estreito em escala log) chega ao ramo final (5/10). O gargalo é o **viés intrínseco do estimador** (S4): para
protegê-lo a camada bilateral custa +0,43 dex de largura, e o IC com recalibração (0,3-0,4 dex) já ocupa metade da separação H_A-H_B (0,80 dex). Um estimador sem viés eliminaria `F` e, com ele, a maior parte desta perda (J8).
(c) **Consequência de desenho a decidir (J18):** com o gatilho 1 como rebaixador, um resultado verdadeiramente Newtoniano (nenhum sinal) vira `NAO_DISTINGUE[GATILHO1]` + relatório de debunker, e não `BOTH_FALSIFIED`; o Estágio 2 o tratava como gatilho de escrutínio, não de rebaixamento.
(d) As colunas de T6 (sem L2/L4) e do estimador "corte" estão em `analysis_summary.log`; o corte é robusto em viés (S6) mas sofre da mesma limitação de largura.

## 6. Gatilhos obrigatórios de checagem adversarial de descoberta de nulo

Além da reexecução adversarial padrão (`AGENTS.md` passo 7, agente separado, implementação do zero, **sempre** obrigatória) e do
papel de debunker (`METHODOLOGY_EXTENSIONS.md` §5, obrigatório para qualquer veredito diferente de `NAO_DISTINGUE`/`NAO_EXECUTAVEL`).
**Efeito dos gatilhos na decisão (§5.2):** os gatilhos 4 e 5 são implementados como *gates* (L1) e rebaixam o veredito; os gatilhos 1, 2 (estrato alto, ou estrato baixo na borda superior), 6, 7, 9 e 10
também rebaixam (L2); os gatilhos 3 e 8 apenas **exigem** escrutínio adicional do debunker. Nenhum gatilho pode promover um veredito.

1. *(herdado, Estágio 2 gatilho 1)* `a0` conjunto fora, por mais de uma ordem de grandeza, de AMBOS `a0_A` e `a0_B`.
2. *(herdado, adaptado, gatilho 2)* Para **cada** estrato: `converged_bracket = False` ou `f̂_k` a menos de `xtol` de `f_lo` ou `f_hi`.
   `f̂_1 = 0` no estrato baixo é **esperado** se esse estrato é limpo (na simulação: 10% das réplicas com `f_baixo=0,02` — cfg 1 — e 60% sem contaminação — cfg 12; T3); a regra
   de fronteira é J5. `f̂_2` na borda, ou `f̂_1` na borda superior, sempre aciona o gatilho.
3. *(herdado, adaptado, gatilho 3)* `f̂_2` fora da faixa `0,25–0,47` da literatura (`PREREGISTRATION_STAGE2.md:35,380`) por margem grande.
   A comparação só vale para o estrato RUWE-alto; para o baixo a faixa não se aplica (J5). Atenção a I-9 (Notas): não verifiquei se a
   faixa de Chae vem de amostra com corte de RUWE.
4. *(substitui o gatilho 4 do Estágio 2)* **Gate de universalidade** (E5) reprovado: `T > T*` ou `zmax > zmax*`. Racional: gravidade (qualquer
   `a0`) prediz o mesmo `δ_b` depois do âncora em todos os estratos; um artefato dependente de RUWE prediz `Δ_b ≠ 0`. É o discriminador
   que separa `M_conv` de gravidade (Seção 3). O gatilho 4 antigo olhava só o padrão *agregado* e, por desenho, não via isto
   (`ADVERSARIAL_DEBUNKER_REPORT.md:249-256`).
5. **Gate de qualidade do ajuste** (E6) reprovado: `χ² > χ²*`. (Retrospectivamente, os números publicados do Estágio 2 dariam
   `χ²≈45,6` para 4 bins no `a0` ajustado — `synthetic_identifiability.log`, linhas `RETRO-GOF` — o que o Estágio 2 não testava.)
6. *(herdado, gatilho 5)* O veredito mecânico depende de `B` ou da seed do bootstrap (refazer com outra seed e `B/2`).
7. **Heterogeneidade ausente.** Se `Δf` e `D` têm IC incluindo 0, o resultado contradiz `DISC-CLAIM-008`: investigar bug/ amostra
   errada antes de qualquer veredito (`DISC-CLAIM-008` registra IC do subgrupo RUWE-alto inteiramente acima de zero).
8. **Composição de RUWE por bin de `g_N`.** Tabela `n_{2,b}/n_b` (b=0..4), só da descoberta. A simulação (configs 10/11) mostra que é neste
   caso — fração RUWE-alta variando com `g_N` — que o `a0` escalar fica enviesado e o condicional o corrige. Se a fração varia pouco,
   o condicional não deve mudar o `a0`; uma mudança grande sem variação de composição aponta para artefato, não para multiplicidade.
9. **Discordância entre estimadores.** `|log10(â0_conjunto / â0_corte)|` maior que a largura combinada dos ICs (ou > 0,3 dex, J6).
10. **Sensibilidade pré-declarada (4 variantes, contadas na Seção 7):** (a) borda de RUWE em 1,3 e 1,5; (b) `K=3` com bordas {1,2; 1,4};
    (c) ruído astrométrico do mock inflado por `max(1, RUWE)` (prática Gaia para RUWE>1,4); (d) `n_mc` e `B` reduzidos à metade.
    Qualquer variante que mude a categoria do veredito ⇒ `NAO_DISTINGUE`.

## 7. Correção para comparações múltiplas (contagem declarada ANTES de olhar qualquer resultado real)

Contagem total declarada: **2 + 3 + 4 = 9** comparações/verificações, de três naturezas diferentes:

| família | quantas | o que são | tratamento |
|---|---|---|---|
| **Confirmatória** | **2** | H_A e H_B testadas contra o **mesmo** IC95% conjunto de `a0` (um ajuste, 2×5 pontos) | sem correção adicional, **igual a** `PREREGISTRATION_STAGE2.md:530-541` e `PREREGISTRATION.md:411-419`; se o autor preferir Bonferroni (IC97,5% para as duas), J11 |
| **Gates** | **3** | heterogeneidade (E2, informativa), universalidade (E5, com 4 bins; omnibus `T` e `zmax` — um só critério), qualidade do ajuste (E6) | cada gate só pode **rebaixar** o veredito para `NAO_DISTINGUE`, nunca promover H_A/H_B; portanto múltiplos gates tornam o procedimento mais conservador, não menos (α nominal 0,05 por gate; limiares empíricos J6) |
| **Sensibilidade** | **4** | variantes (a)-(d) do gatilho 10 | descritivas; qualquer mudança de categoria rebaixa a `NAO_DISTINGUE`; nunca usadas para escolher a variante favorável |

Sem busca sobre número de bins de `g_N`, sobre o bin-âncora (=4), nem sobre os valores de `f_k` além das K bisseções. **Não há** seleção entre
estimadores depois do resultado: o estimador primário é o conjunto (J2); o corte e os ajustes por estrato são descritivos.
O número de execuções sobre dado real autorizadas por este rascunho, se um dia travado, é **uma** (mais as 4 variantes declaradas);
qualquer execução adicional exige novo pré-registro.

## 8. Linha vermelha: o holdout selado NÃO é tocado neste Estágio, sob nenhuma circunstância

Esta seção herda `PREREGISTRATION_STAGE2.md:441-470` e `METHODOLOGY_EXTENSIONS.md` §6 e os **endurece** em dois pontos
(8.3 e 8.5), onde o registro anterior deixava uma folga.

**8.1 Regra.** Os **12.944** sistemas do holdout (`discovery_holdout_split.json`, seed=20260814; `PREREGISTRATION.md:100`)
permanecem **intocados** pelo Estágio 3: nem lidos, nem usados para calibrar `f_k`, ajustar `a0`, calibrar limiares de gate,
estimar contagens por estrato, ou "checar sanidade". Isto vale para este rascunho, para o seu eventual LOCK, e para a
execução sobre a amostra de descoberta que o LOCK venha a autorizar.

**8.2 Uso futuro exige lock formal próprio e autorização humana separada.** Nenhum item deste documento — nem o LOCK de
`PREREGISTRATION_STAGE3.md`, nem a conclusão de uma análise sobre a descoberta, nem uma reexecução adversarial bem-sucedida —
equivale a autorização para abrir o holdout. Qualquer abertura futura exige, **cumulativamente**:
(i) uma entrada própria no `DECISION_LEDGER.yaml` que trave o holdout para uma análise nomeada (hipótese, estatística,
regra de decisão e código congelados por hash antes da abertura); (ii) autorização humana explícita e separada do autor;
(iii) passagem pelo Gate de Replicação (`03_REPLICATION_GATE/PROTOCOL.md`, onde o holdout é "finalmente aberto e testado —
nunca antes"; `METHODOLOGY_EXTENSIONS.md` §6). A frase de `TEST_QUEUE.yaml:409-411` ("holdout ... permanece intocado,
**disponível** para essa tentativa futura") **não** deve ser lida como autorização — ver inconsistência I-3 em `STAGE3_DRAFT_NOTES.md`.

**8.3 Regras de código (mais estritas que as do Estágio 2).**
- Nenhum script do Estágio 3 pode acessar a chave `holdout_pair_ids` nem `n_holdout` de `discovery_holdout_split.json`
  (nem por `json.load` seguido de indexação, nem por iteração sobre `.items()`/`.keys()` do dicionário inteiro — o Estágio 2
  só proibiu a indexação literal, o que não impede varredura do dicionário). Carregar o JSON deve ser feito por função que
  devolve **apenas** `discovery_pair_ids` e descarta o dicionário imediatamente.
- O grep de portão (Estágio 2, Seção 9.3) deve ser estendido: além de `catalog.parquet` e das chaves de holdout, procurar
  `.items()`, `.keys()`, `.values()`, `list(split)` e `**split` aplicados ao dicionário do split.
- **Contagens, quantis, médias e qualquer agregado** só podem ser calculados **depois** de filtrar pela lista de descoberta;
  nenhum agregado sobre a tabela completa (por exemplo "fração de RUWE>1,4 no arquivo") pode ser calculado ou logado.

**8.4 Calibração de limiares de gate e de contagens mínimas** usa **somente** simulação sintética (Seção 9) e a amostra de
descoberta; nunca o holdout. O tamanho documentado do holdout (12.944) é usado neste rascunho **apenas como um número de N**
em simulações sintéticas (configs 5, 8 e 13) para saber qual poder um dia teria uma análise com esse N — isto não lê nenhuma
chave do holdout.

**8.5 Lacuna a fechar antes de qualquer LOCK (J14, Seção 11).** Pelo que está documentado, `quality_filtered_sample.parquet`
contém as 43.147 linhas (descoberta + holdout) e o Estágio 2 o carregou inteiro antes de filtrar. Para que a "linha vermelha"
signifique "nenhuma linha do holdout em memória", e não só "nenhuma coluna do holdout computada", propõe-se (decisão do autor)
um **extrato só-descoberta** criado em passo separado, com sha256, e usado como única entrada do Estágio 3. Enquanto isso não
for decidido, a garantia do Estágio 3 é a mesma, mais fraca, do Estágio 2 — e deve ser descrita assim no relatório de resultados.

**8.6 Se um passo do Estágio 3 parecer exigir o holdout** (por exemplo, "para ter mais N no estrato RUWE-alto"): o passo é
**cancelado**, o veredito do Estágio 3 é `NAO_EXECUTAVEL`/`NAO_DISTINGUE` (Seção 5) e a necessidade é registrada como argumento
para uma decisão humana separada. Nunca "só olhar".

## 9. Checklist de revalidação sintética obrigatória imediatamente antes do dado real

Mesma lógica de `PREREGISTRATION_STAGE2.md:472-528` e `PREREGISTRATION.md` §4b/4c/5b: revalidar sob os **parâmetros exatos** que vão
tocar o dado, ANTES de tocá-lo. Todos os itens abaixo usam **somente** simulação sintética.

- **9.1 Integridade dos módulos.** `sha256` de `selfcal_pipeline.py`, `companion_injection.py`, `build_synthetic_population.py`,
  `deprojection_common.py` registrado no commit de lock; `fit_a0()` conserva a guarda `a0>0` (`PREREGISTRATION_STAGE2.md:483-486`).
- **9.2 Reexecutar `synthetic_identifiability.py` sem alteração de lógica** com os valores finais de `n_mc`, `B`, `xtol` e com
  **≥ 200 réplicas nulas** (configs 0, 1 e 12) para recalibrar `T*`, `zmax*`, `χ²*` (J6) e o fator `F` (J7). Critérios C1-C3 da Seção 5.1
  devem continuar passando.
- **9.3 Testes de equivalência U1-U2** (gerador condicional == `build_fake_real_dataset` no caso escalar; pipeline com `f` vetorial
  constante == escalar) passam com diferença exatamente 0.
- **9.4 Validações A e B do Estágio 1 (escalar) reexecutadas** sob os novos `n_mc`/`B` (`PREREGISTRATION_STAGE2.md:488-496`) e **análogos
  condicionais**: população homogênea ⇒ `Δf` não detectado (≤10%); população heterogênea ⇒ `f̂_k` recuperados dentro de tolerância
  pré-declarada.
- **9.5 Stress de mis-especificação** (configs 3, 4 e 13): taxa com que os gates detectam o artefato e viés do `a0` conjunto vs. corte,
  nos valores finais; se o gate de universalidade tem poder < 0,5 ao N real, **documentar** (Seção 5.3) — não corrigir depois.
- **9.6 Checagem de insensibilidade do bin-âncora ao `a0`**, **por estrato** (`_check_anchor_bin_mond_negligible`,
  `validate_b_recover_a0_with_contamination.py:59-66`).
- **9.7 Ausência de dado real e do holdout** por grep estendido (Seção 8.3) em **todos** os `.py` do Estágio 3 (incluindo os do protótipo).
- **9.8 Registro de seeds** e do `sha256` de cada script e log; seeds do Estágio 3 real fora do bloco reservado.
- **9.9 Revisão de código do FPB por agente separado** (o ponto de maior risco é o descasamento entre as duas execuções pareadas por estrato —
  mesma classe de risco de `PREREGISTRATION_STAGE2.md:284-294`): confirmar que `raw_a` e `raw_b` usam a mesma seed e a mesma
  ordem de sistemas.
- **9.10 Portão de contagem (E0) como primeiro passo real**, com a regra de parada da Seção 4.2.

**Critério de bloqueio:** falha em qualquer item ⇒ o Estágio 3 não toca dado real; a falha é documentada antes de qualquer nova tentativa.

## 10. O que NÃO está sendo testado

- **Tamesis.** Nenhum resultado deste Estágio — positivo, negativo ou inconclusivo — confirma, refuta, detecta ou favorece
  Tamesis, nem o favorece sobre ΛCDM ou MOND (`AGENTS.md`, Proibições). `a0_A = cH0/2π` reproduz uma coincidência numérica já
  conhecida na literatura MOND (Milgrom, arXiv:2001.09729, citado em `PREREGISTRATION.md:78-84`): H_A sobreviver não é evidência
  específica de Tamesis. A nota de identificabilidade de `PREREGISTRATION_STAGE2.md:76-83` aplica-se sem diluição.
- **Nenhum Problema do Milênio** é resolvido, aproximado ou tornado alcançável por este teste ou por este rascunho.
- **Que MOND/AQUAL seja a gravidade correta.** Só se testa se o `a0` ajustado sob a função `ν` "simple" fixa é compatível com
  `a0_A` ou `a0_B`, depois de uma correção de contaminação.
- **A causa física do RUWE alto** (companheira oculta genuína vs. solução astrométrica degradada vs. cruzamento de ambos).
  `DISC-CLAIM-008` já registra essa lacuna (`CLAIM_LEDGER.yaml`, `known_gaps`). O modelo `f_multi(RUWE)` é uma **parametrização
  estatística de contaminação**, não uma afirmação de que RUWE alto = companheira.
- **Dependência de `f_multi` em outras covariáveis** (massa, separação/`g_N`, distância, cor, magnitude) nem dependência
  contínua de `f_multi` em RUWE **dentro** de um bin. A Seção 0 mostra que a primeira é exatamente o caso em que o modelo
  escalar falha e o condicional ajuda — e é desconhecida.
- **Os herdados do Estágio 1:** `γ_M = -0,7` fixo, wobble de fotocentro simplificado, função `ν` "simple", desprojeção de Chae.
- **Qualquer dado do holdout** (Seção 8) e **qualquer replicação independente** (Gate de Replicação).
- **Um modelo de mistura hierárquico/MCMC** (a "Linha G" do `DIAGNOSTICO_AUTOMACAO_LINHAS_PESQUISA.md:60-66`). Seria uma frente
  distinta, com pré-registro próprio (nota datada em `DIAGNOSTICO_AUTOMACAO_LINHAS_PESQUISA.md:128-141`).
- **Um corte de exclusão "RUWE<1,4" como método de medição** — ele aparece aqui apenas como estimador **alternativo** de
  comparação (J2). O texto do diagnóstico ("fixar ... o corte de RUWE < 1.4") descreve outra frente.

## 11. Julgamentos ABERTOS para o autor humano resolver ANTES de qualquer LOCK (NÃO resolvidos por este rascunho)

Cada item traz as opções, o que a simulação (Seção 0) informa e **o que ela não informa**. A ordem é por consequência.

**J1. Papel da amostra de descoberta: exploratória ou confirmatória?** A forma `f_multi(RUWE)`, a borda 1,4 e o estrato "RUWE-alto" foram motivados
por uma checagem do debunker **sobre a própria amostra de descoberta** (`ADVERSARIAL_DEBUNKER_REPORT.md:304-335`). Opções: (a) tratar uma execução na descoberta como
**exploratória** (veredito limitado a `NAO_DISTINGUE`/descrição) e reservar o confirmatório a um holdout sob lock próprio (Seção 8.2); (b) tratá-la como confirmatória
(não recomendado: dupla utilização); (c) não executar na descoberta. A aritmética retrospectiva (Seção 0) indica que os gates provavelmente reprovariam na descoberta, o que
torna (a)/(c) mais naturais. *A simulação não resolve J1: é uma questão de desenho inferencial.*

**J2. Estimador primário.** (a) ajuste conjunto condicional; (b) apenas o estrato baixo ("corte" RUWE≤1,4); (c) apenas diagnóstico (nenhum veredito sobre H_A/H_B).
A simulação: nos cenários consistentes com o modelo (configs 0, 1) os três dão o mesmo viés; com artefato aditivo não modelado no estrato alto (configs 3, 4) o conjunto e o escalar
pioram e o corte não; com composição variando com `g_N` (configs 10, 11) o escalar piora e o condicional/corte não. O corte custa ~19% da amostra e **não** elimina a
contaminação residual do estrato baixo (`f̂_1` não-nulo). *Não informado:* se o estrato baixo real é limpo.

**J3. Bins de RUWE.** K=2 {1,4} (proposto); K=3 {1,2; 1,4}; forma contínua; RUWE por componente em vez de `max(RUWE1,RUWE2)`. K=3 exige contagens por célula que só o portão E0 revela;
não testei K=3 sinteticamente (apenas declarado como variante de sensibilidade, gatilho 10b). A distribuição real de RUWE por bin de `g_N` é desconhecida para mim (V2).

**J4. Contagem mínima por bin (`N_A_MIN`, `N_S_MIN`) e regra se violada.** A simulação sugere `N_A_MIN ≈ 500` sistemas **no bin-âncora do menor estrato** (Seção 5.1) e a mesma ordem de grandeza
em cada célula `(k,b)`; escolha final e a regra (parar vs. fusão pré-declarada de estratos) são do autor. Se a fração RUWE-alta varia com `g_N` (configs 10/11), as células de
canto ficam muito menores do que `N·π/5`.

**J5. Fronteira de `f̂_k`.** `f̂_1 = 0` (estrato baixo limpo) acontece com frequência (Seção 0, T3) e é fisicamente esperado (`f ≥ 0`); a regra do pipeline travado devolve a extremidade de menor |δ| e marca
`converged_bracket=False`. Declarar: (a) fronteira em 0 no estrato baixo é permitida (só sinalizada) e `f̂_2` na fronteira nunca; (b) como o FPB trata réplicas que cruzam 0 (aqui: clip em 0).
O resíduo do âncora do estrato baixo **a f=0** reportado pelo debunker é **negativo** (−0,012; `ADVERSARIAL_DEBUNKER_REPORT.md:323`): se confirmado, a calibração do estrato baixo não terá troca de sinal.

**J6. Gates: forma, α, limiares.** Omnibus `T` vs. `zmax` com Bonferroni vs. χ²; limiar nominal (χ²) vs. empírico. Na simulação o bootstrap do procedimento é **conservador**: o desvio-padrão de `Δ_b` no FPB excede a variação entre réplicas por fator tipicamente ~2 (faixa 1,2-5,2 nas configs; T5), então os `z` são subestimados e
gates com limiar nominal (χ²) teriam falso positivo abaixo do nominal e menos poder; os limiares empíricos aqui vêm de apenas 30 réplicas nulas. Recalibrar com ≥200 réplicas (Seção 9.2).

**J7. Camada de viés.** Estágio 2: dividir o IC por 1,4–1,6 (unilateral). Proposta: fator bilateral `F3` derivado da simulação (Seção 5). Ou substituir por calibração por simulação (J8).

**J8. Viés intrínseco do estimador de `a0` (+~0,2 dex).** Mesmo sem contaminação (config 12) e em populações homogêneas, o `a0` recuperado excede o injetado por fator ~1,6 e o IC95% (largura ~0,18 dex
no N da descoberta) exclui o valor injetado em todas as réplicas. Opções: (a) aceitar e usar a camada de viés (com a regra completa, H_A verdadeiro nunca chega a `A_ONLY` — 0/120, T8 — e H_B verdadeiro só em 5/10; o método só consegue, na prática, falsificar ou declarar `NAO_DISTINGUE`);
(b) corrigir o estimador antes de qualquer lock (p.ex. ajustar `a0` comparando `δ` observado com `δ` **simulado pelo próprio pipeline** em vez de `δ_AQUAL` analítico em `g_N` projetado);
(c) não travar. *Não investiguei a causa raiz* — apenas medi o viés. Ver I-8 em `STAGE3_DRAFT_NOTES.md` (o Estágio 1 atribuiu o viés ao wobble; o controle positivo de SPARC-004 já mostrava razão 1,44 sem contaminação).

**J9. Parâmetros de execução real:** `n_mc` (200?), `B` (2000?), `xtol` (5×10⁻⁴), **seeds do Estágio 3 real** (fora do bloco 20260955000-999) e a convenção de sub-seeds derivadas do pipeline (D3 nas Notas);
propagação da incerteza de `f̂` no IC (FPB com recalibração, como aqui, ou `f̂` fixo, como no Estágio 2 `bootstrap_a0_refit`). Na simulação o IC com recalibração é ~2× mais largo
que com `f̂` fixo (T2) — escolha afeta diretamente a cobertura.

**J10. Pesos do ajuste conjunto:** contagens (aqui), inversa da variância bootstrap, ou iguais (Estágio 2). Não testei alternativas.

**J11. Nível do IC / Bonferroni para H_A e H_B** (95% como no Estágio 2, ou 97,5% para 2 hipóteses).

**J12. Como implementar `f` vetorial.** O protótipo faz *monkeypatch* em tempo de execução de `ci.sample_companion_injection_mc` (nenhum arquivo editado). Para o dado real: parâmetro aditivo em
`selfcal_pipeline.py` (como `return_raw`, DISC-DEC-029) — **altera módulo travado**, exige decisão — ou código novo em arquivo novo que chama o pipeline por estrato (não precisa de `f` vetorial: o
conjunto só usa `δ` por estrato).

**J13. Variante de sensibilidade com ruído astrométrico do mock inflado por RUWE** (gatilho 10c): incluir ou não; conta como comparação (Seção 7).

**J14. Extrato "só-descoberta"** (Seção 8.5): criar (passo de preparação com autorização própria) ou aceitar a garantia mais fraca do Estágio 2.

**J15. Estratificador adicional/alternativo** (`astrometric_excess_noise`, `ipd_*`): exigiria dado novo (Seção 2.2).

**J16. Abrir `catalog.parquet` bruto** para o baseline populacional de RUWE pré-corte: exigiria autorização separada (Seção 2.2).

**J17. Escopo do Estágio 3:** (a) teste H_A/H_B com gates; (b) diagnóstico do estrato RUWE-alto, sem veredito sobre H_A/H_B; (c) arquivar sem lock. A recomendação do rascunho está na Seção 0.

**J18. Zona de "não distingue":** o rascunho define-a mecanicamente (Seção 5.3); decidir se a largura máxima do IC (hoje: `log10(a0_B/a0_A)` = 0,80 dex) e os gatilhos 1, 2, 6, 7, 9, 10 são aceitos como **rebaixadores** automáticos. Consequências medidas (T8): (i) o ramo de largura (L4) elimina quase todo `A_ONLY` com H_A verdadeiro porque o alargamento bilateral absorve o viés intrínseco (J8); (ii) o gatilho 1 rebaixa um resultado Newtoniano a `NAO_DISTINGUE[GATILHO1]` em vez de `BOTH_FALSIFIED` (7/10 em cfg 2) — no Estágio 2 esse gatilho só exigia escrutínio. O autor pode preferir tornar o gatilho 1 apenas de escrutínio.

**J19. Escolha dos botões sintéticos** (`f_low=0,02`, `f_high=0,45`, jitter 0,3 km/s, composição 0,04→0,41): foram escolhidos por mim para ser plausíveis/ilustrativos, e o jitter após uma exploração (parte `explore`).
O autor pode exigir outros cenários antes do lock (checklist 9.2/9.5).

## 12. Proveniência da simulação sintética (tudo que sustenta os números da Seção 0)

| arquivo (em `stage3_ruwe_conditional/`) | função |
|---|---|
| `synthetic_identifiability.py` | script final (v2: v1 + modo `--part merge` + lista padrão de configs 0-13) |
| `synthetic_identifiability_v1_used_in_run1.py` | cópia **exata** do script usado nas execuções run1 e run2 |
| `synthetic_identifiability_run1_cfg0-9.log/.json` | execução 1: partes `unit`, `explore`, `retro`, e `main` para configs 0-9 (100 réplicas), 18:04-18:56 UTC |
| `synthetic_identifiability_run2_cfg10-13.log/.json` | execução 2: `main` para configs 10-13 (40 réplicas), 19:07-19:28 UTC |
| `synthetic_identifiability.log` / `.json` | **fusão** (`--part merge`) de run1+run2: cópia do texto das partes unit/explore/retro, todas as 140 linhas de réplica com seeds, resumos por config, gates, características operacionais |
| `analyze_synthetic_identifiability.py` → `analysis_summary.log` | tabelas T1-T8 (viés, largura/cobertura, detecção, testes de duas amostras, gates, características operacionais T6, regra completa T8, varredura em N T7) lidas **só** do JSON fundido |

```text
a05f5af188e178e244b30f7c7420631181c5e790e5d79402115309d2db5c0438  synthetic_identifiability.py
3bfe8ef9b9bf57e13843d4e100a41cc7b92ee2bb81cbf5f7742b68efab9ff450  synthetic_identifiability_v1_used_in_run1.py
7cec26333874acf5c8ee5eb5b67c064b30d1f50c8cd21cc8c67778753dbfce72  analyze_synthetic_identifiability.py
6e65ebfc231935f69c2b6594e29756063ffc45735191a251b4b4f04a307458e5  synthetic_identifiability_run1_cfg0-9.json
444db5b13271ac4077181d44a05ea57c14bc952451a002a3a027dfaf91e96652  synthetic_identifiability_run1_cfg0-9.log
b02b4aeed8af3f58f2f6f6bc35c1ef8cd1d78b99a4902b792a686e7cb6895ccf  synthetic_identifiability_run2_cfg10-13.json
ba9c1ede92417f69cbeb0eaffb3a8d2002bfb2f984065e72f2528509062a8e43  synthetic_identifiability_run2_cfg10-13.log
0c8fcf8266f27c71b551abbc2afaa8464542de8c420738394ea07b2c7e5887b8  synthetic_identifiability.json
3c624c8c473395db775948b695e818483ba04cdef22e342a87e5d9a6c78a243a  synthetic_identifiability.log
61363602cdc60e42017914d55258f2a0ac52c8a94cca21651085037e55974f32  analysis_summary.log
```

**Por que duas execuções + fusão (registro honesto, Notas D4).** A lista padrão de configs da execução 1 só continha 0-9 (erro meu: as configs 10-13 foram adicionadas ao dicionário `CONFIGS` mas não ao padrão de `--configs`). Percebi quando a execução terminou com 100 réplicas.
Em vez de refazer as 100 réplicas, rodei só as configs 10-13 com o **mesmo código** (v1) e acrescentei ao script o modo de fusão (somente adição de código; nenhuma lógica de simulação alterada — a diferença v1→v2 é exatamente a função `part_merge`, a opção `--merge`/`choices` e o padrão `--configs`).
Como o código é determinista por seed, rodar tudo de uma vez daria o mesmo resultado; não verifiquei isso rodando (custaria ~100 min).

**Ambiente:** Python 3.11.17, numpy 2.4.6, scipy 1.17.1 (venv no scratchpad; a máquina não tinha numpy). Comando: `python -I -B synthetic_identifiability.py --part {all|main|merge} ...` com `sys.path` explícito e `sys.dont_write_bytecode=True` (nenhum `__pycache__` criado; conferido).

**Parâmetros:** `N_MC=100`, `B=300`, `xtol=2×10⁻³`, passo de `f` do FPB 0,15, `f∈[0;0,9]`, fração do estrato alto `5784/30203`, `a0_A=1,082288×10⁻¹⁰`, `a0_B=6,800218×10⁻¹⁰`, bins e âncora como no Estágio 2. 14 configs × 10 réplicas. Para o dado real: `N_MC=200`, `B=2000`, `xtol=5×10⁻⁴` (J9).

**Seeds (registro):** réplica = `20260955000 + 6·(10·config + rep) + papel`, papéis 0-5 = catálogo, verdade, calibração escalar, calibração estrato baixo, calibração estrato alto, bootstrap. Usadas: **840 valores distintos, de 20260955000 a 20260955839** (listadas no campo `seeds_used_main` do JSON e na linha de cada réplica do log).
Dev/unit: 20260955960-20260955967 (testes U1-U5) e 20260955968-20260955973 (parte `explore`, reusadas em todos os botões). Nenhuma seed fora do bloco 20260955000-999 foi passada por este script.
**Ressalva:** o pipeline travado deriva sub-seeds por deslocamentos fixos (`+111_111_111`, `+222_222_222`, `+333_333_333`, `+444_444_444`, `+555_555_555`, `+777`), que não estão no bloco numérico literal (Notas D3).

**Testes unitários (log, parte `unit`):** U1 `build_truth(f escalar)` ≡ `build_fake_real_dataset`: diferença máxima **0** (PASS); U2 pipeline com `f` vetorial constante ≡ escalar: **0** (PASS); U3 erro da interpolação linear em `f` do FPB: ≤ **0,0039 dex**;
U4 população homogênea `f=0,30`, N=6000: `f̂=0,3249`, convergido; U5 contagens iguais por bin: `[1200]×5`. (O catálogo sintético **não** equalizado dá `[2910, 1019, 690, 559, 817]` em N=6000.)

**Exploração do botão de artefato (log, parte `explore`; 1 réplica, N=12.944, seeds 968-973):** `Δ_b` (bins 0-3) por botão —
escala do wobble ×4/×8/×16 dá `Δ` **negativo** no bin 0 (`-0,24/-0,45/-0,67`: o `f̂` do estrato alto absorve e supercorrige); escala do ruído ×3/×5/×8 dá `Δ_bin0 = -0,00/+0,04/+0,10` e `Δ_bin1 = +0,07/+0,11/+0,17`, com bins 2-3 ≤ 0,04 (mais fraco que a assinatura do debunker);
o jitter aditivo de velocidade dá `Δ>0` decrescente com `g_N` (0,2 km/s: `[0,12; 0,21; 0,05; 0,04]`; 0,3: `[0,24; 0,29; 0,10; 0,06]`; 0,5: `[0,37; 0,39; 0,17; 0,12]`), a assinatura do debunker — por isso 0,3 km/s foi escolhido para as configs 3 e 13 e 0,1 km/s para a config 4.

**Como reproduzir (sem tocar dado real):**
```text
python -I -B synthetic_identifiability_v1_used_in_run1.py --part all --workers 4 --configs 0,1,2,3,4,5,6,7,8,9 --outdir <dir> --prefix run1
python -I -B synthetic_identifiability_v1_used_in_run1.py --part main --workers 4 --configs 10,11,12,13 --outdir <dir> --prefix run2
python -I -B synthetic_identifiability.py --part merge --merge <dir>/run1.json <dir>/run2.json --outdir <dir> --prefix synthetic_identifiability
python -I -B analyze_synthetic_identifiability.py --json <dir>/synthetic_identifiability.json --out <dir>/analysis_summary.log
```

---

## [Preenchido depois da análise] Resultado

*(vazio — este documento é um rascunho; nenhuma análise sobre dado real foi executada. Só pode ser preenchido após LOCK em ação humana separada e execução autorizada em separado.)*

## [Preenchido depois da reexecução adversarial] Veredito adversarial

*(vazio — mesma razão.)*

---

## Notas e correções datadas do referee hostil (2026-10-25, `DISC-DEC-156`)

Fonte: `referee/REFEREE_REPORT.md` (reprodução independente de 4 células: 134 campos cada, diferença 0 — prova determinismo/integridade do JSON, NÃO a correção do estimador). Verificação da sessão: a aritmética de F1 foi refeita (0,798 − 2·log10(1,65) = 0,363 dex < 2b ≈ 0,42 dex para b ≈ +0,21 dex).

- **Correção (F1, alta):** a regra de decisão T8 não consegue emitir `A_ONLY` nem `BOTH_FALSIFIED` (incompatibilidade algébrica entre `lo ≤ a0_A` e o ramo L4 quando o viés é ≥ ~0,2 dex). "0 falsificações falsas em 130 réplicas" é vácuo. A regra é um teste **sem poder** para H_A, H_B e o caso Newtoniano; o rótulo "Segurança" não se sustenta.
- **Correção (F2, alta):** o texto (§8.4, E0, item 9.10) abre brecha para dado real/leitura do split inteiro antes de qualquer LOCK, contradizendo o §8.1 e `PREREGISTRATION_STAGE2.md:460-470`. Default correto sob J14 indeciso: `NAO_EXECUTAVEL`. O código da simulação está limpo (sem acesso a holdout/parquet).
- **Correção (F3):** `N_A_MIN ≈ 450` vem do critério C3 (largura < 0,80 dex), incompatível com L4 (< 0,36 dex); nenhum N testado dá poder sob a regra completa.
- **Nota (F4):** limiares dos gates vêm de 30 nulos (IC90 largos); falso positivo ≈20% em N=12.944 fora da calibração; o erro astrométrico subestimado (mecanismo do debunker) não está na grade principal e passa pelos gates.
- **Nota (F5):** S2 (p=0,88) é tautologia do gerador; a "Razão 2" do §0.4 e a leitura de I-4 não se sustentam nele.
- **Nota (F6–F8):** "o condicional corrige" vale só no ponto (IC médio 2,56 dex em cfg 11); interpolação do bootstrap validada só no regime interior; números decisivos são provisórios e o ramo L6 fica sem limite pré-declarado.
- **Reordenação do §0.4 (substitui a ordem original):** (1) regra sem poder por construção; (2) brechas de holdout no texto; (3) viés e F circulares; (4) 19 julgamentos abertos. Poder "10/10" e "~8/10 em N=12.944" vêm de 10 réplicas in-model e não devem ser citados como evidência de identificabilidade real.
- **Veredito:** o rascunho NÃO é candidato a LOCK. A recomendação "não travar" permanece. Status continua `DRAFT`.
