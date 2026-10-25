# REFEREE_REPORT — árbitro hostil independente sobre `PREREGISTRATION_STAGE3_DRAFT.md` (SPARC-FMULTI-STAGE3) e sua simulação sintética

**Autoridade:** `DISC-DEC-155` (`00_GOVERNANCE/DECISION_LEDGER.yaml:10074-10119`), critérios de rejeição (a)-(d).
**Data:** 2026-10-25. **Escopo respeitado:** nenhum comando git; nenhum dado real baixado/aberto; holdout não tocado; nenhum arquivo existente editado;
scripts do árbitro apenas em `stage3_ruwe_conditional/referee/`; Python do venv do scratchpad (`python -I -B`, sem bytecode; `find` não achou `__pycache__` nem `.whl`).
Seeds usadas pelo árbitro: as do rascunho (cfg 0 rep 0, cfg 10 rep 0, cfg 11 rep 0, cfg 12 rep 3) e, só no teste de linearidade, 20260955900-902 (bloco reservado, não usadas pelo rascunho).

---

## 0. Veredito em uma página

| critério DISC-DEC-155 | veredito do árbitro | resumo |
|---|---|---|
| (a) caminho para tocar holdout/dado real sem lock posterior | **ATINGIDO no texto; NÃO atingido no código** | Os scripts estão limpos (grep + imports). Mas o **texto** do rascunho tem três brechas: §8.4 autoriza "calibrar limiares de gate e contagens mínimas" com "a amostra de descoberta" (dado real, antes de qualquer LOCK); §4.3 E0 / §9.10 inserem um passo de dado real numa seção declarada "somente sintética" e com a frase ambígua "antes do LOCK dos resultados"; e o E0 como escrito exige `read_parquet` da tabela inteira (43.147 linhas, holdout incluído), o que contradiz §8.1 ("nem lidos") e **não** cumpre `PREREGISTRATION_STAGE2.md:460-470` ("nunca materializar a lista de holdout em memória") — ao contrário do que o §8 diz ("endurece"). |
| (b) zona de veredito ambígua/reinterpretável pós-hoc | **ATINGIDO parcialmente** | A zona `NAO_DISTINGUE` é mecânica e exaustiva (bom), mas: todos os números que decidem (`F`, `T*`, `zmax*`, `χ²*`, `N_A_MIN`, `n_mc`, `B`) são provisórios; o ramo L6 é "SEMPRE sujeito a debunker" sem limites pré-declarados (precedente: no Estágio 2 o debunker derrubou o veredito mecânico com uma análise nova pós-hoc); o gatilho 9 tem "ou" ambíguo; L4 e a regra completa foram escritos **depois** de ver as larguras (Notas D6). |
| (c) validação sintética incapaz de distinguir condicional de escalar | **ATINGIDO para a regra ponta-a-ponta; NÃO atingido para o sub-enunciado "detectar heterogeneidade"** | A regra completa nunca emite `A_ONLY` nem `BOTH_FALSIFIED` em 140 réplicas (nem no cenário Newtoniano) e isso é **estrutural** (§2 abaixo), não amostral: é um teste sem poder para H_A. O método distingue heterogeneidade (in-model) e o `a0` condicional difere do escalar só quando a composição de RUWE varia com `g_N` (cenário assumido, e com IC inútil no cfg 11). |
| (d) inconsistência com `PREREGISTRATION_STAGE2.md` / `TEST_QUEUE.yaml` / `CLAIM_LEDGER.yaml` | **NÃO atingido nos fatos e números; ATINGIDO em itens menores de procedimento/vocabulário** | Todos os números e citações que conferi batem. Inconsistências: `:460-470` do Estágio 2 (materialização), vocabulário novo (`NAO_DISTINGUE`/`NAO_EXECUTAVEL`/`A_ONLY`) sem mapeamento para o vocabulário do registro (`INCONCLUSIVO`, `preregistered_inconclusive`), gatilho 1 passa de "escrutínio" a "rebaixador" (declarado, J18). |

**Reprodução:** 4 células reproduzidas **bit a bit** (diferença exatamente 0 em 134 campos cada). Isso prova determinismo e integridade do JSON, **não** correção do estimador (usei o mesmo script; não reimplementei o pipeline).

**Resultado global:** o rascunho, como está, **deve ser rejeitado como candidato a LOCK** (pelos critérios a, b, c); isso coincide com a recomendação do próprio §0.4, que o árbitro considera **correta no sentido, mas fundamentada em parte por argumentos fracos e subestimando o achado mais forte** (§5 abaixo). Como registro de diagnóstico (status DRAFT, sem LOCK), é utilizável **depois** de corrigidos os achados F1-F4.

---

## 1. Reprodução independente (tarefa 1)

Script: `referee/repro_cell.py` (importa `synthetic_identifiability.run_replicate` e roda `(cfg, rep)` com as seeds registradas por `seed_for`, `n_mc=100`, `n_boot=300`);
comparação: `referee/compare_repro.py` → `referee/compare_repro.out`.

| célula | seeds (cat..boot) | campos comparados | diferenças ≠ 0 | exemplos (registrado = reproduzido) |
|---|---|---|---|---|
| cfg 0 rep 0 (homogêneo, H_A, N=30203) | 20260955000-005 | 134 | **0** | `a0` escalar 1,759002170468684e-10; `a0` cond. 1,7850581103604303e-10; IC recal [1,2517836e-10; 2,5489073e-10]; `T`=1,96559541 |
| cfg 10 rep 0 (composição cresce com `g_N`) | 20260955600-605 | 134 | **0** | `a0` escalar 5,671273934737432e-11; cond. 1,869888129837024e-10; `f̂`(baixo, alto) = (0,018849; 0,430752) |
| cfg 11 rep 0 (composição decresce com `g_N`) | 20260955660-665 | 134 | **0** | `a0` escalar 3,9220380439032063e-10; cond. 1,5253271841555806e-10; IC recal [2,99e-13; 6,82e-10] |
| cfg 12 rep 3 (sem contaminação) | 20260955738-743 | 134 | **0** | `a0` escalar 1,8087381403012383e-10; cond. 1,7641822212424937e-10 |

Também conferi: sha256 de todos os arquivos listados no §12 do rascunho **batem** (10/10); `diff` v1→v2 do script é só `part_merge` + opção `--merge` + padrão de `--configs` (como declarado); 840 seeds distintas de 20260955000 a 20260955839; testes U1/U2 = 0 exatamente, U3 = 0,00392 dex, U4 `f̂`=0,32488, U5 = [1200]×5 (batem com o log). A tabela §5.5 do rascunho (regra Estágio 2 e regra Estágio 3 por config) foi **recomputada do JSON por mim** (`referee/check_claims.py`) e bate integralmente.

Observação: o rascunho (§12) diz "não verifiquei rodando que rodar tudo de uma vez daria o mesmo resultado". Para as 4 células amostradas, rodar célula-a-célula dá o mesmo resultado que o JSON fundido; isso sustenta (não prova para as 140) o argumento de determinismo.

## 2. Afirmações numéricas verificadas (tarefa 2)

### 2.1 Contra JSON/log (tudo conferido; nenhuma diferença numérica material)
S1 (detecção 10/10, 10/10, 10/10, **2/10** em N=3000; FP 0/10 em cfg 0, 8, 12; `D` 10/10, 10/10, 10/10, 8/10) ✔. S2 (+0,216±0,055 vs +0,213±0,028; Welch p=0,882; cfg 12 +0,212±0,035) ✔. S3 (cfg 10: −0,252 / +0,229; cfg 11: +0,595 / +0,195±0,108 / corte +0,230; Estágio 2 `BOTH_FALSIFIED` 20/20) ✔. S4 (larguras, cobertura 0/10; 30/30 `[VIES]` no Estágio 2 em cfg 0, 1, 12) ✔. S5 (0,18 / 0,33-0,40 / 0,18-0,20; 0,62; 1,02; 2,10) ✔. S6 (cfg 3: +0,350/+0,409/+0,201, p=9,8e-6 e 2,0e-7 e 0,46; cfg 13: +0,351/+0,386/+0,205) ✔. S7 (limiares `T*`=2,4516, `zmax*`=1,1119, `χ²*`=5,2037; cfg 3 10/10/10; cfg 13 7/8/8; cfg 4 1/2/4) ✔. S8 (`F3`=1,65339; 130 réplicas = 13 configs com `a0` definido; tabela §5.5) ✔ — **mas ver F1**. S9 (`Δ_b`=[0,229; 0,219; 0,187; 0,102], `z`=[4,6; 4,7; 4,2; 2,3], `T`=66,15; `χ²`=45,55 no `a0_fit` publicado; a tabela do debunker `:319-323` e `RESULTS_PRIMARY_STAGE2.md:95-99` conferem com os insumos) ✔.
`5784/30203 = 0,191504` ✔ (`ADVERSARIAL_DEBUNKER_REPORT.md:179`; `RESULTS_PRIMARY_STAGE2.md:225-226`: 19,15%).

### 2.2 Contra os arquivos citados (amostra de ~30 citações `arquivo:linha`; todas bateram)
`TEST_QUEUE.yaml:403-411` (`next_action`: reabertura condicional a RUWE, "holdout ... disponível") ✔; `PREREGISTRATION_STAGE2.md:35,64-70,87,124,134-138,303-307,326,333-353,380,441-445,460-470` ✔; `RESULTS_PRIMARY_STAGE2.md:25,32,45-50,78,95-99,192-200,225-228,246-249` ✔; debunker `:94-100,122,164-173,179,192-195,249-256,304-335,319-323,350-358,470-474` ✔; `PREREGISTRATION.md:8-18,65-72,78-84,94-102,411-419` ✔; `run_stage2_primary_analysis.py:105,130-145,177-205` ✔. **Não** conferi: `PREREGISTRATION_STAGE2.md:179-222,250-294,472-541`, `METHODOLOGY_*`, `PROVENANCE_CHAE_EQS.md`, `DIAGNOSTICO_*` (fora do que o árbitro leu).

### 2.3 As inconsistências I-1..I-11 das Notas (o rascunho lista "10 inconsistências com o registro")
I-1 ✔ (`PREREGISTRATION_STAGE2.md:1,3,549,678`). I-2 ✔ (`RESULTS_PRIMARY_STAGE2.md:45-50` diz "Dois dos cinco" e lista um; tabela `:192-200` e `CLAIM_LEDGER` dizem "só o gatilho 3"). I-3 ✔ (citação literal de `TEST_QUEUE.yaml:409-411`). **I-4: o fato é correto, mas a inferência é exagerada** (ver F5). I-5 ✔ **e mais forte do que o rascunho admite**: `run_stage2_primary_analysis.py:178` lê o parquet inteiro, `:181-185` carrega o JSON de split inteiro, o `merge` com Hwang exige `len(merged)==len(qf)` e só depois filtra por `isin(discovery_ids)` — se `qf` tivesse só 30.203 linhas o filtro seria desnecessário; logo `qf` tem as linhas do holdout em memória. I-6, I-7, I-9, I-10: plausíveis, não verifiquei na fonte. I-8: o fato empírico (cfg 12, sem companheiras: +0,212 dex, indistinguível de cfg 0/1, p=0,96/0,86) ✔; ressalva: o "1,44" do controle positivo (`PREREGISTRATION.md:11-12`) é razão de **sinal δ no bin 0**, o 1,63 é razão de **`a0`** — grandezas diferentes (esclarecimento F11). I-11 ✔ (nenhum resíduo `.whl` encontrado).

### 2.4 Regra de decisão T8 — 130 réplicas, 0 falsificações falsas, 0/120 `A_ONLY` (o ponto central)

`referee/check_claims.py` → `referee/check_claims.out`. Recomputei a regra completa do JSON: coincide com a coluna "Estágio 3 proposto" da §5.5. Totais sobre as 140 réplicas: `ND[VIES]` 30, `ND[GATE]` 42, `ND[LARG]` 43, `ND[GATILHO1]` 7, `ND[CI]` 10, `ND[GATILHO2]` 3, **`B_ONLY` 5, `A_ONLY` 0, `BOTH_FALSIFIED` 0.**

**A regra é incapaz de chegar a `A_ONLY` — e por razão algébrica, não por amostragem.**
- `A_ONLY` exige `lo ≤ a0_A` (S_A). Com viés do estimador `b ≈ +0,21 dex` e IC aproximadamente simétrico em log, isso exige meia-largura inferior ≥ `b`, i.e. largura `w ≳ 2b ≈ 0,42 dex`.
- O ramo L4 exige `W' = w + 2·log10 F < log10(a0_B/a0_A)` = 0,7982, i.e. `w < 0,7982 − 0,4367 = 0,3614 dex`.
- Como `F = 10^max|viés|` foi **calibrado nesse mesmo viés** (log10 F = 0,2184), as duas condições só são compatíveis se `b < 0,7982/4 = 0,1995 dex`. O viés observado (cfg 0, 1, 12: +0,2184/+0,2154/+0,2124 medianas; réplicas individuais nas 30 nulas: +0,10 a +0,27; o mínimo individual entre as 120 réplicas H_A é +0,05, em cfg 11, com IC de 3 dex) é maior na média; a conclusão estrutural vale para o IC típico, e o resultado observado confirma: nenhuma réplica satisfaz ambas as condições. **Observado: 0/120 réplicas H_A satisfazem `lo ≤ a0_A` e `w < 0,361` simultaneamente** (cfg 0: lo≤a0_A 0/10, w<0,36 10/10; cfg 1: 5/10 e 1/10; cfg 5: 8/10 e 0/10; …; nunca ambos).
- Logo "0/120 `A_ONLY`" **não é um achado empírico com 10 réplicas por config; é uma identidade do desenho** (`F` bilateral + L4 + viés ≈ 0,21). O rascunho (§5.5(b), J18) atribui a frequência à "largura do IC alargado" e sugere que "um estimador sem viés eliminaria `F`" — correto, mas **não diz que, com o estimador atual, a regra é matematicamente incapaz de dar suporte a H_A**.

**"0 falsificações falsas em 130" é vácuo.** A regra proposta **nunca emitiu `BOTH_FALSIFIED` em nenhuma das 140 réplicas**, incluindo as 10 do cenário Newtoniano (cfg 2), onde `BOTH_FALSIFIED` é a resposta correta e o Estágio 2 acertava 10/10; a regra proposta deu `ND[GATILHO1]` 7/10, `ND[GATE]` 2/10, `ND[LARG]` 1/10. O conjunto de saídas alcançáveis é {`NAO_DISTINGUE[*]`, `B_ONLY`}. Uma regra que nunca falsifica nada trivialmente "nunca falsifica falsamente". O limite de Wilson 95% superior para 0/130 é 2,9%, mas isso é irrelevante: a taxa de falsificação **verdadeira** também é 0.
- A palavra "Segurança" em §5.5(a) e "ao custo de poder total" em S8 são honestas apenas pela metade: o texto precisa dizer com todas as letras que **esta regra, aplicada a qualquer dado com o estimador atual, não pode dar suporte a H_A nem falsificar H_A, H_B ou o caso Newtoniano; ela só pode devolver `NAO_DISTINGUE` ou (assimetricamente) `B_ONLY`** (5/10 em cfg 9, porque o viés **positivo** empurra o IC para `a0_B`). A assimetria contra H_A (a hipótese ligada à linha Tamesis) é uma propriedade do desenho que o registro deveria declarar, mesmo que não seja intencional.
- Contrafactual (**in-sample, circular, só ilustrativo**; `referee/counterfactual_onesided.py`): com a mesma máquina de gates/gatilhos mas correção de viés **unilateral** (`IC/F`, como o Estágio 2) em vez de bilateral+L4, saem 77 `A_ONLY` (de 120 H_A), 10/10 `B_ONLY` em cfg 9, 1 `BOTH_FALSIFIED`. Ou seja: **todo o "poder" ou "falta de poder" é decidido pela convenção de tratar o viés, e `F` vem do próprio gerador (inverse crime)** — qualquer das duas escolhas converte um viés não explicado de um gerador sintético em veredito sobre dado real. Isto reforça que J7/J8 é pré-requisito lógico, não detalhe.

## 3. Circularidade e desenho pós-hoc — quanto invalida cada conclusão (tarefa 3)

| conclusão | circularidade/ajuste | quantificação do árbitro | efeito |
|---|---|---|---|
| S1 poder de detectar heterogeneidade | gerador e estimador compartilham o modelo (Chae Eqs. 11-13, mesmo wobble, mesmo ruído) — "inverse crime", Notas A1/R1; só 3 mis-especificações, escolhidas pelo autor | `D` real do debunker = +0,237 (`:335`); cfg 1 dá +0,276 com `f_low=0,02`, `f_high=0,45` — **dois botões ajustados para acertar um número**; "reproduz D≈+0,28" (§0.4 item 2) não é teste | S1 vale como **limite superior** do poder real. Não invalidada, mas não transportável; além disso a heterogeneidade já é observada no dado (IC do debunker exclui 0), então S1 tem pouco valor decisório |
| S2 (`a0` escalar não enviesado por heterogeneidade) | **tautologia do gerador**: com rótulo RUWE independente de tudo, a população agrupada é exatamente iid Bernoulli(`f̄`); os p-valores (0,88) testam uma nulidade verdadeira por construção. IC95% da diferença cfg1−cfg0 ≈ ±0,041 dex | zero informação além de "o código implementa o que diz" | **não** serve de evidência para a realidade; depende de A2 (RUWE ⟂ `g_N`, massa, distância), que o próprio S9 (Δ_b≠0, z≈4,6) contradiz |
| S3 (condicional corrige quando a composição varia com `g_N`) | cenário H2 é **assumido** (π por bin) — não há evidência de que exista no dado (V2) | cfg 11: largura do IC recal. **2,56 dex em média; mediana 3,28; 6/10 réplicas >3 dex** (cobertura 10/10 trivial); n_âncora(alto)=236. Só cfg 10 tem IC informativo (0,29 dex) e ainda assim o veredito é `ND[VIES]` 8/10 | corrige o **ponto**, que ainda carrega +0,2 dex; não corrige inferência |
| S4 viés intrínseco +0,21 dex / `F`=1,65 | `F=10^max|viés mediano|` calculado em cfg 0, 1, 12 (as mesmas 30 réplicas que calibram os gates); causa raiz desconhecida (Notas R7, J8); pode ser específica do gerador | 13/13 configs com `a0` definido têm viés positivo no condicional (mín. +0,162, cfg 9); 43% das 30 réplicas nulas individuais têm viés > log10 F | `F` não é uma constante do método; é uma constante deste gerador |
| S7 gates | limiares = percentil 95 de **30** nulos; taxas nesses nulos são in-sample | IC90 bootstrap do limiar: `T*` [1,89; 6,79], `zmax*` [0,93; 1,85], `χ²*` [2,73; 6,43]. FP "qualquer gate" in-sample 3/30 (Wilson [0,03; 0,26]); fora-da-amostra (bootstrap cruzado) média 0,12, IC90 [0; 0,30]; **leave-one-config-out: 0,20 em cada config nula**. Fora do conjunto de calibração: **cfg 8 (homogêneo, N=12.944) dispara 3/10, cfg 5 (N=12.944) 1/10 → FP ≈ 4/20 = 20% em N=12.944** com limiares calibrados em N=30.203 (`referee/gate_sensitivity.py`). Poder cfg 4 (0,1 km/s): 0,5 (IC90 [0,1; 0,6]); cfg 13: 0,9 [0,8; 0,9] | poder dos gates "10/10" vale para o artefato **escolhido por exploração para parecer com o debunker** (0,3 km/s) |
| S7 / gates vs. mecanismo físico | a **mesma** exploração (`explore`, 1 réplica, N=12.944, no JSON/log) mostra que subestimar o ruído astrométrico do estrato alto (o mecanismo que o debunker cita, `:164-173`) dá viés condicional **+0,31 (×3) e +0,34 (×5) dex com os três gates silenciosos** (`T`=0,3/0,8 < 2,45; `χ²` 1,3/2,0 < 5,2); essa família foi descartada da grade principal (cfgs 3, 4, 13 só têm o jitter aditivo) e não aparece em §0.3 S7 | 1 réplica por botão: indicativo, não estatístico, mas é a omissão mais perigosa do desenho de mis-especificação |
| Limiares/L2/L4/regra completa | escritos depois de ver T6/T2 (Notas D6, D2); catálogo de contagens iguais, aparo do bin 0 (mediana ajustada ao valor real publicado) e botão 0,3 km/s foram escolhidos depois de execuções de desenvolvimento (que existem no scratchpad: `dev_explore*.out`) | não quantificável; direção: tudo o que foi ajustado melhora a **verossimilhança** do cenário, nada melhora o poder da regra | o desenho final **não** é cego às execuções; o rascunho declara isso (D2, D6) — OK em honestidade, mas implica que nenhuma taxa de §5.5 é uma estimativa fora-da-amostra |
| FPB com recalibração linear em `f` (§4.3 E3) | "erro ≤ 0,004 dex" vem de U3: f 0,10→0,25, N=1500, uma realização, estrato de baixo `f` | `referee/linearity_check.py` (estrato alto, f_true=0,45, N≈2556, par 0,40/0,55): erro linear máx. **0,018 dex em t=+0,33**, 0,067 (t=−0,67), 0,165 (t=−1,33), 0,277 (t=−2), 0,247 (t=+2), **0,602 (t=+3,3)**. O FPB extrapola até t≈3 quando `f` recalibrado cruza o intervalo (cfg 11: sd(f̂ alto)=0,049; n_âncora=236) | larguras de IC (S5, L4) e `sd_FPB` (que é ~2× a variação entre réplicas, T5, sem explicação) dependem de uma aproximação validada só no regime interior; a largura de 2,56-3,3 dex do cfg 11 pode ser parcialmente artefato numérico (não demonstrado) |

## 4. Disciplina do holdout (tarefa 4)

**Código (critério a, nível script): limpo.**
- `grep -n -i "holdout_pair|catalog\.parquet|\.parquet|split\.json|discovery_pair|pandas|data/|quality_filtered|hwang"` nos três `.py` do Estágio 3 → **nenhuma ocorrência**. "holdout" aparece só em docstring (l.28), na constante `N_HOLDOUT_DOCUMENTED = 12944` (l.86) e em nomes de config (`..._N12944`).
- Imports: `deprojection_common`, `companion_injection`, `selfcal_pipeline`, `build_synthetic_population`, `delta_obs_newt` (só `BIN_EDGES_LOG_GN_SPARC003`); nenhum importa pandas ou abre parquet/JSON de dado (`grep` de `read_parquet|open(` nos módulos: só docstrings). `analyze_synthetic_identifiability.py` só lê o JSON sintético.
- O uso de `12944` como N de simulação não lê nenhuma chave do holdout. A premissa A3 (`π_alto` do holdout = o da descoberta) é uma suposição, não um dado.

**Texto do rascunho (critério a, nível documento): três brechas** (F2 abaixo). Ordem de gravidade:
1. §8.4: "Calibração de limiares de gate e de contagens mínimas usa **somente** simulação sintética (Seção 9) **e a amostra de descoberta**". Isso autoriza ajustar limiares com dado real **antes** de qualquer LOCK e contradiz §0, §5.2 (limiares "fixados antes do LOCK"), §5.3(i). É a única frase do rascunho que permite explicitamente tocar dado real fora do passo E0.
2. §4.3 E0 ("primeiro e único acesso a dado real **antes do LOCK dos resultados**") e §9.10 (passo E0 como item do checklist "**todos os itens abaixo usam somente simulação sintética**"): o texto não diz se o E0 ocorre antes ou depois do LOCK do pré-registro; as Notas §7(a) dizem "depois do LOCK", o rascunho não. E0 é um passo de dado real dentro de uma seção declarada sintética.
3. O E0 como escrito ("contagens da AMOSTRA DE DESCOBERTA já filtrada pelos 30.203 ids") exige carregar a tabela completa e o dicionário de split inteiro (I-5). §8.3 propõe "função que devolve apenas `discovery_pair_ids` e descarta o dicionário imediatamente" — mas `json.load` **materializa** a lista de holdout em memória, e `PREREGISTRATION_STAGE2.md:467-470` exige "nunca materializar a lista de holdout em memória nem em log". O §8 diz que "endurece" o Estágio 2; neste ponto ele **não cumpre** o Estágio 2. §8.5/J14 reconhece a lacuna, mas o default (J14 não decidido) é "a mesma garantia mais fraca" — o default deveria ser `NAO_EXECUTAVEL`. Além disso, criar o "extrato só-descoberta" também lê a tabela completa uma vez (J14 não resolve a materialização, só a move para um passo separado).
4. Menor: o grep de portão proposto (§8.3) não impede `json.load(...)` seguido de `del`; só varredura de dicionário. Precisa de parse em streaming ou de extração por processo separado que escreva só `discovery_pair_ids`.

O texto fecha corretamente: `catalog.parquet` bruto (§2.2/J16, autorização separada), uso futuro do holdout (§8.2: ledger próprio + autorização humana + Gate de Replicação), e a leitura de `TEST_QUEUE.yaml:409-411` ("disponível") como **não** autorização (I-3). Esses pontos estão bem.

## 5. A recomendação "não travar" (§0.4) está justificada? (tarefa 5)

**A conclusão (não travar) está justificada — por mais razões do que o rascunho dá, e as razões listadas têm pesos desiguais.**

*O que a evidência sustenta com firmeza:*
- Com o estimador atual (`b≈+0,21`), a regra proposta **não pode** devolver suporte a H_A nem falsificação (§2.4) — conclusão **estrutural**, mais forte que o "0/120" do rascunho. Isso sozinho inviabiliza o LOCK como "teste de H_A/H_B".
- Há lacunas de texto sobre o holdout (§4) e 19 julgamentos abertos, vários estruturais; números decisivos são provisórios. LOCK é impossível sem eles.
- O viés de +0,21 dex sem nenhuma companheira (cfg 12) é um achado robusto no gerador e relevante para o registro (contradiz a causa atribuída no Estágio 2 §7 — "wobble"); a causa raiz permanece desconhecida e pode ser do gerador. Isso é **mais forte** do que a Razão 2 que o rascunho enfatiza.
- Execução na descoberta ⇒ `NAO_DISTINGUE[GATE]` quase certo (S9): `z=4,6` vs `zmax*`≈1,1 (mesmo com sd inflado 2× e limiar duplicado, ainda dispara). Aritmética grosseira (ignora covariâncias e usa `f̂` escalar para os dois estratos), mas a margem é grande.

*O que a evidência **não** sustenta (ou sustenta menos do que o texto sugere):*
- Razão 2 ("a premissa motivadora não se traduz em viés de `a0`", S2): depende de uma tautologia do gerador e de A2 (RUWE independente de tudo), que S9 contradiz. Não deveria figurar como razão de peso.
- Qualquer número de **poder** ("10/10", "~8/10 em N=12.944", "n_âncora ≳ 230") com 10 réplicas, limiares calibrados com 30 nulos em N=30.203 e FP de gate ≈20% em N=12.944: são ordens de grandeza in-model, não estimativas. O limiar "N_A_MIN≈450" vem de C3 (`w<0,80`), inconsistente com L4 (`w<0,36`); nenhum N testado satisfaz L4 + `lo≤a0_A` (F3).
- A frase "viés intrínseco do **método**" (S4): é do estimador **no gerador** (inverse crime); sem causa raiz (J8) não se sabe se existe no dado real.
- A Razão 4 ("com N do holdout o IC mal cabe em 0,80 dex") usa IC com recalibração linear não validada em extrapolação (§3).

Em resumo: **diagnóstico "não travar" — correto; classificação de gravidade do rascunho — subestima F1 (estrutural) e F2 (texto do holdout), superestima S2/Razão 2 e os números de poder.** O árbitro recomenda reordenar §0.4: (1) regra sem poder por construção; (2) brechas do holdout; (3) viés/F circular; (4) julgamentos abertos.

## 6. Achados classificados e ranqueados (tarefa 6)

Legenda: **ERRO** = requer correção no texto/desenho; **ESCL** = esclarecimento/nota. Severidade: ALTA / MÉDIA / BAIXA.

| # | classe | sev. | critério | achado | correção/nota proposta |
|---|---|---|---|---|---|
| F1 | **ERRO** | **ALTA** | (c) | A regra completa (§5.2) **não pode** devolver `A_ONLY` nem `BOTH_FALSIFIED` com o estimador atual: `F` bilateral calibrado no viés + L4 ⇒ condição `b<0,1995 dex` não satisfeita (b≈0,21). 0/120 `A_ONLY`, 0/140 `BOTH_FALSIFIED` (inclusive Newtoniano), 5 `B_ONLY` só com H_B verdadeiro (assimetria). "0 falsificações falsas em 130" é vácuo; §5.5(a) chama isso de "Segurança" | Escrever em §0.3/S8, §5.5 e §5.3: "esta regra é um teste **sem poder** para H_A, H_B e Newtoniano; saídas alcançáveis: `NAO_DISTINGUE` e (assimetricamente) `B_ONLY`"; remover "Segurança"; declarar que o resultado é algébrico; não usar a regra como candidata a LOCK antes de J7/J8 |
| F2 | **ERRO** | **ALTA** | (a) | Brechas de texto sobre dado real/holdout: §8.4 ("e a amostra de descoberta" para calibrar limiares); E0/§9.10 (passo de dado real em seção "somente sintética"; "antes do LOCK dos resultados" ambíguo); E0 exige carregar tabela/JSON inteiros (holdout em memória) em desacordo com §8.1 e com `PREREGISTRATION_STAGE2.md:460-470`; default de J14 = garantia fraca | Remover "e a amostra de descoberta" de §8.4; mover E0/9.10 para "após LOCK + autorização"; E0 só com extrato só-descoberta (sha256) produzido por processo separado com parse em streaming; default sem extrato = `NAO_EXECUTAVEL`; corrigir a afirmação "endurece" |
| F3 | **ERRO** | MÉDIA-ALTA | (c),(b) | `N_A_MIN≈450` / C3 (`w<0,80`) incompatível com L4 (`w<0,36`): nenhum N testado dá poder sob a regra; C1/C3 "passam" por estimativa pontual de 10 réplicas (FP 0/10 tem Wilson sup. 0,28; poder 8/10 inf. ~0,49) | Derivar `N_A_MIN` do critério efetivo da regra (alcançabilidade de `A_ONLY` e `B_ONLY` sob L4) ou remover; reportar IC de Wilson em C1-C3 |
| F4 | **ERRO** (omissão) | MÉDIA | (c) | S7 "os gates detectam esse artefato": só para o artefato escolhido para ser detectável. Na própria exploração (JSON) o artefato "erro astrométrico subestimado" (mecanismo do debunker) dá viés +0,31/+0,34 dex com gates silenciosos. FP de gate fora do conjunto de calibração: cfg 8 3/10, cfg 5 1/10 (N=12.944); LOO 20%; IC90 dos limiares muito largo | Acrescentar a S7 e §5.5: poder dos gates por **família** de artefato (incluindo ruído); declarar limiares como não transportáveis entre N e parâmetros (`n_mc=200`, `B=2000`); manter exigência ≥200 nulos **por N de interesse** |
| F5 | ESCL | MÉDIA | (c) | S2 é tautologia do gerador (iid Bernoulli de `f̄`); depende de A2, contradita por S9. I-4 e Razão 2 de §0.4 se apoiam nisso | Reclassificar S2 como "verificação de código/identidade", retirar p-valores como evidência e a Razão 2 como argumento de peso; manter I-4 como ponto de lógica a levar ao autor, não como inconsistência do registro |
| F6 | ESCL | MÉDIA | (c) | S3 "o condicional o corrige": no cfg 11, IC médio 2,56 dex (mediana 3,28), cobertura 10/10 trivial; ponto ainda +0,2 dex. Só cfg 10 tem IC informativo e o veredito é `ND[VIES]` | Qualificar S3: "corrige o ponto, não a inferência" |
| F7 | ESCL | MÉDIA | (b),(c) | Interpolação linear em `f` do FPB validada só no interior (U3, 0,004 dex); no estrato alto o erro vai de 0,018 (t=0,33) a 0,6 dex (t≈3,3) e o FPB extrapola até t≈3; `sd_FPB` ≈2× a variação entre réplicas, sem explicação | Corrigir "≤0,004 dex" em §4.3 E3 (condicionar ao regime U3); item novo no checklist §9: validar linearidade no estrato alto e limitar `|t|` ou recalibrar por bisseção; explicar o fator ~2 antes de usar `sd_FPB` em gates |
| F8 | ESCL | MÉDIA | (b) | Zona de veredito: números decisivos provisórios; L6 "SEMPRE sujeito ao debunker" sem limites; gatilho 9 com "ou" ambíguo; ordem de avaliação dos ramos só implícita; L4/regra completa escritos depois de ver larguras (D6) | Declarar que o debunker só **rebaixa** e só com checagens pré-declaradas (novas análises exigem novo pré-registro — precedente do Estágio 2); fixar gatilho 9 numa só forma; escrever "primeiro ramo que dispara, na ordem L0..L6"; listar "números a fixar antes do LOCK" como bloqueio |
| F9 | ESCL | MÉDIA | (d) | `D` do cfg 1 "reproduz" o do debunker porque `f_low`/`f_high` foram escolhidos para isso (2 botões, 1 alvo); diferença 0,276 vs 0,237 (+16%) | Reescrever "reproduz" como "foi calibrado para se aproximar de" |
| F10 | ESCL | BAIXA-MÉDIA | (d) | Vocabulário novo (`NAO_DISTINGUE`, `NAO_EXECUTAVEL`, `A_ONLY`, `B_ONLY`, `BOTH_SURVIVE`) sem mapeamento para o registro (`INCONCLUSIVO`, `H_A_FALSIFICADA`, `evidence_level: preregistered_inconclusive`); gatilho 1 passa de escrutínio (Estágio 2 `:357`) a rebaixador (declarado J18) | Tabela de mapeamento no rascunho; manter J18 como bloqueio de LOCK |
| F11 | ESCL | BAIXA | (d) | I-8 compara "1,44" (razão de sinal δ no bin 0) com "1,63" (razão de `a0`); o ponto de fundo (cfg 12) é válido | Ajustar a redação de I-8 |
| F12 | ESCL | BAIXA | — | Incidente de escopo D1 (`pip download` criou um `.whl` fora do diretório permitido; removido). Conferido: nenhum resíduo. DISC-DEC-155 diz "Arquivos criados apenas em stage3_ruwe_conditional/" — violação transitória, já declarada | Manter registro D1; sem ação |
| F13 | ESCL | BAIXA | — | Seeds: sub-seeds do pipeline (`+111_111_111`, …) fora do bloco literal (I-7/D3) — inevitável sem editar módulo travado | Manter ressalva |
| F14 | ESCL | BAIXA | — | R11/V7: DISC-DEC-155 prevê "spot-check independente da sessão antes do referee"; as Notas dizem que nenhuma verificação independente ocorreu até a escrita. Este relatório não substitui esse spot-check de procedência | Registrar se o spot-check foi feito pela sessão |

## 7. O que este relatório não verificou

- Não reimplementei o estimador/pipeline (a reprodução usa o mesmo código); não investiguei a causa raiz do viés +0,21 dex.
- Não rodei as 140 réplicas (4 células; ~minutos cada, 0 diferenças). O teste de linearidade é uma realização (uma seed).
- Não abri `quality_filtered_sample.parquet`, `catalog.parquet` nem o split (proibido); a conclusão sobre as 43.147 linhas em memória é inferência a partir do código de `run_stage2_primary_analysis.py` e da documentação, como no rascunho.
- Citações `arquivo:linha` conferidas por amostragem (~30); `METHODOLOGY_EXTENSIONS.md`, `PROVENANCE_CHAE_EQS.md`, `03_REPLICATION_GATE/PROTOCOL.md`, `DIAGNOSTICO_*` e `PREREGISTRATION_STAGE2.md:179-222,250-294,472-541` não foram relidos.

## 8. Arquivos do árbitro (todos em `referee/`)

| arquivo | função |
|---|---|
| `REFEREE_REPORT.md` | este relatório |
| `repro_cell.py` | reprodução de uma célula `(cfg, rep)` com as seeds registradas |
| `compare_repro.py` / `compare_repro.out` | diff campo a campo (134 campos × 4 células; 0 diferenças) |
| `check_claims.py` / `check_claims.out` | recomputo da regra T8, alcançabilidade de `A_ONLY`, limiares LOO, contagens do Estágio 2 |
| `gate_sensitivity.py` / `gate_sensitivity.out` | sensibilidade do poder/FP dos gates à incerteza do limiar |
| `counterfactual_onesided.py` / `counterfactual_onesided.out` | contrafactual in-sample com correção de viés unilateral (ilustrativo, circular) |
| `linearity_check.py` / `linearity_check.out` | erro da interpolação/extrapolação linear em `f` no estrato alto |
