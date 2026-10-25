# Notas do rascunho de pré-registro `SPARC-FMULTI-STAGE3` (`DISC-DEC-155`)

**Status:** notas de trabalho do agente de rascunho; acompanham `PREREGISTRATION_STAGE3_DRAFT.md` (Status: DRAFT).
**Data:** 2026-10-25. **Nenhum dado real, nenhum holdout, nenhum comando git.** Escrito só em `stage3_ruwe_conditional/`.
Conteúdo: o que li e fiz (§1), o que assumi (§2), o que não consegui verificar (§3), desvios e incidentes (§4), inconsistências com o registro (§5), riscos (§6),
mapa para o critério de refutação de `DISC-DEC-155` (§7) e o que o humano decide primeiro (§8).

## 1. O que foi lido e feito (rastro)

**Leitura integral (antes de escrever qualquer linha):** `DECISION_LEDGER.yaml` (DISC-DEC-155, linhas 10074-10119; também DEC-023/026/027/029/032),
`PREREGISTRATION_TEMPLATE.md`, `AGENTS.md`, `METHODOLOGY_EXTENSIONS.md`, `PREREGISTRATION.md` (SPARC-004), `PREREGISTRATION_STAGE2.md`
(e o cabeçalho de `PREREGISTRATION_STAGE2_DRAFT.md`), `RESULTS_PRIMARY_STAGE2.md`, `RESULTS_SUMMARY_STAGE1.md`, `METHODOLOGY_ADDENDUM.md`,
`PROVENANCE_CHAE_EQS.md`, `adversarial_stage2/ADVERSARIAL_DEBUNKER_REPORT.md`, `analysis/selfcal_pipeline.py`, `build_synthetic_population.py`,
`companion_injection.py`, `deprojection_common.py` (importado), `validate_b_recover_a0_with_contamination.py` + `results/validation_B_run.log`,
trechos de `TEST_QUEUE.yaml` (SPARC-004), `CLAIM_LEDGER.yaml` (DISC-CLAIM-008), `DIAGNOSTICO_AUTOMACAO_LINHAS_PESQUISA.md` (Linha G),
`03_REPLICATION_GATE/PROTOCOL.md` (início) e `run_stage2_primary_analysis.py` (só as funções de carga e o grep de portão, para a Seção 8).

**Não abri:** nenhum `.parquet`, `.json` de dado ou de split, nenhum arquivo de `data/` (apenas listei o diretório `COSMOLOGY_WIDE_BINARIES_MC_DEPROJECTION/`),
nenhum resultado `results/*.json` de dado real. **Nenhum comando git.** Nenhum download de dado.

**Simulação:** `synthetic_identifiability.py` (sha256 do script e dos logs na Seção 12 do rascunho),
Python 3.11.17 / numpy 2.4.6 / scipy 1.17.1 em venv no scratchpad (a máquina não tinha numpy). Execução com `python -I -B`, `sys.path` explícito,
`sys.dont_write_bytecode = True` (nenhum `__pycache__` fora do diretório permitido; verificado por `ls`).

## 2. O que assumi

A1. **A "verdade" sintética usa o mesmo modelo de contaminação que o pipeline** (injeção de Chae Eqs. 11-13 + wobble simplificado). Logo os testes de
    recuperação são favoráveis ("inverse crime"); as únicas mis-especificações testadas são as 3 que eu escolhi (escala do wobble, escala do ruído, jitter aditivo de velocidade).
A2. **RUWE entra como um rótulo binário de estrato** (≤1,4 / >1,4), sorteado independentemente de tudo, exceto nas configs 10/11 (onde a fração do estrato alto depende do bin de `g_N`).
    Não modelei RUWE contínuo, nem o RUWE por componente, nem correlação com massa/distância.
A3. **Fração do estrato RUWE-alto = 5784/30203** (`ADVERSARIAL_DEBUNKER_REPORT.md:179`), assumida igual para amostras de qualquer N (inclusive o N do holdout).
A4. **`f_low=0,02`, `f_high=0,45`** são botões escolhidos para reproduzir qualitativamente o padrão do debunker (estrato alto com `f` na faixa 0,25-0,47, estrato baixo quase limpo, prevalência
    agrupada ≈0,10≈`f̂` do Estágio 2). **Não são estimativas do dado real.**
A5. **Catálogo estrutural:** gerador original do Estágio 1 (separação log-uniforme etc.), com duas adaptações declaradas: (i) contagens iguais por bin de `g_N` (a amostra real tem ~6.040 por bin,
    `PREREGISTRATION_STAGE2.md:127`; o gerador original dá 2910/1019/690/559/817 em N=6000 — U4); (ii) cauda inferior do bin 0 aparada até a mediana igualar a do real já publicada
    (1,8807×10⁻¹⁰, `RESULTS_PRIMARY_STAGE2.md:95`; sem o ajuste a mediana sintética era ~5,7×10⁻¹¹ e o sinal MOND no bin 0 ficava superestimado). As medianas dos bins 1-4 sintéticos
    ficaram dentro de ~22% das reais (a pior é o bin-âncora: 5,9×10⁻⁸ vs 4,86×10⁻⁸).
A6. **Parâmetros de custo:** `n_mc=100` (real: 200), `B=300` (real: 2000), `xtol=2×10⁻³` (real: 5×10⁻⁴), 10 réplicas por config. Os ICs simulados são, portanto, mais ruidosos que os reais,
    e as taxas estimadas com 10 réplicas têm IC de Wilson largo (reportados).
A7. **Recalibração no bootstrap por interpolação linear em `f`** entre duas execuções pareadas (erro medido em U3: ≤0,004 dex num caso de teste), em vez de refazer a bisseção por réplica.
A8. **Limiares dos gates** (`T*`, `zmax*`, `χ²*`) = percentil 95 das 30 réplicas nulas modelo-consistentes (configs 0, 1, 12) — grosseiro; as taxas nessas configs são in-sample.
A9. **Pesos do ajuste conjunto = contagens por (estrato, bin)** e fator de viés bilateral `F3 = 10^max|viés mediano|` nas configs modelo-consistentes — escolhas minhas, marcadas como julgamentos abertos.
A10. **Estimador alternativo "corte"** = estrato baixo isolado, método escalar travado (sem camada de gates).

## 3. O que NÃO consegui verificar

V1. Se `quality_filtered_sample.parquet` contém de fato as 43.147 linhas (descoberta + holdout) e as colunas `RUWE1`/`RUWE2` — só a documentação, não o arquivo (proibido abrir).
V2. A distribuição real de RUWE por bin de `g_N` (se a fração do estrato alto varia com `g_N`; é exatamente o que decide entre "condicional = escalar" e "condicional ≠ escalar", configs 10/11 vs 1).
V3. Se a faixa `f_multi = 0,25-0,47` e os valores 0,48/0,36 de Chae vêm de amostras com corte de RUWE (o debunker cita "69,2% com RUWE<1,2", `ADVERSARIAL_DEBUNKER_REPORT.md:192-195`, de forma ambígua). Não refiz
    a leitura de Chae (não era o escopo) e não inventei referência alguma: tudo que cito de literatura está citado **por intermédio** dos arquivos do repositório.
V4. Se o excesso do estrato RUWE-alto do dado real é (i) mais multiplicidade do que o modelo captura, (ii) artefato astrométrico, ou (iii) mistura. A simulação só mostra que um artefato aditivo de velocidade
    (configs 3/4) **pode** reproduzir qualitativamente o padrão — não que seja a causa.
V5. A aritmética retrospectiva (`RETRO-GOF`, `RETRO-UNIVERSALIDADE`) usa σ = (IC95% largo)/3,92 sobre tabelas publicadas (N_boot=500 do debunker, N_boot=2000 do Estágio 2), ignora covariâncias e a incerteza de `f̂`.
    É ordem de grandeza, não um teste.
V6. Que o patch em tempo de execução de `ci.sample_companion_injection_mc` (f vetorial) seja o modo certo de implementar o estimador no dado real; é um protótipo (J12).
V7. Nenhuma verificação independente (segundo agente) do protótipo foi feita — só os testes unitários U1-U5. O `STAGE3_DRAFT_NOTES.md` e o rascunho **não** passaram por referee hostil.

## 4. Desvios de procedimento / incidentes (para registro honesto)

D1. **Escrita fora do diretório permitido (transitória):** ao testar disponibilidade de pip, executei `pip download numpy` com o diretório de trabalho em `COSMOLOGY_WIDE_BINARIES_MC_DEPROJECTION/`;
    isso criou um `.whl` (~16 MB) ali. Removi imediatamente (`rm`, conferido por `ls`), mas o mtime do diretório mudou (17:34). Nenhum outro arquivo existente foi tocado.
D2. **Exploração antes da execução final:** rodei execuções de desenvolvimento (N pequeno ou 1-4 réplicas; um teste de fumaça de N=3.000; duas rodadas de 4 réplicas em N=30.203; duas varreduras de botões de 8 combinações em N=12.944) para depurar o código e escolher o botão do artefato;
    isso alterou o desenho (adicionei, **depois** de ver resultados de desenvolvimento: o catálogo de contagens iguais, o aparo do bin 0, as configs de composição por `g_N`, a de artefato aditivo, a config sem contaminação e os gates). Nenhum número de desenvolvimento é citado; a escolha do botão está
    documentada no log final (parte `explore`, seeds 968-973) e as saídas de desenvolvimento ficaram no scratchpad (fora do repositório, não salvas como deliverable). As execuções de desenvolvimento usaram seeds do bloco (as próprias seeds das réplicas finais ou 960-973); nenhuma fora do bloco.
    Isto é flexibilidade de análise legítima num estudo de simulação exploratório, mas deve ser conhecida: o desenho final foi ajustado depois de ver execuções preliminares.
D3. **Seeds:** todas as seeds passadas pelo script estão no bloco 20260955000-999 (registradas no JSON e no log). O pipeline travado deriva sub-seeds por deslocamentos fixos (`+111_111_111`, ..., `+555_555_555`, `+777`), que caem **fora** do intervalo numérico literal — são funções deterministas das
    seeds registradas, mas uma leitura literal de "somente o bloco" não é cumprível sem editar o pipeline travado (proibido).
D4. **Duas coisas aconteceram com as execuções finais:** (i) a primeira execução foi abortada ~1 min depois de iniciada e reiniciada do zero porque três comentários de linha citada no script estavam errados (corrigi só os comentários; nenhum resultado foi visto da execução abortada);
    (ii) a execução reiniciada ("run1") cobriu só as configs 0-9 porque esqueci de acrescentar 10-13 ao valor padrão de `--configs` — percebi ao terminar. Rodei então as configs 10-13 ("run2") com o mesmo código e acrescentei ao script o modo `--part merge`
    para fundir os dois JSONs (somente adição de código; v1 preservado como `synthetic_identifiability_v1_used_in_run1.py`). Detalhes e hashes na Seção 12 do rascunho.
D5. **`pkill -f` acidental:** ao abortar a primeira execução, usei `pkill -f synthetic_identifiability.py` num comando cujo próprio texto continha esse nome, o que encerrou o shell da ferramenta (sem efeito nos arquivos). Registrado por completude.
D6. **A regra completa (T8) foi implementada no pós-processamento depois de ver T6:** a §5.2 (ramos L2 e L4 inclusive) foi redigida depois de eu ver T6 (que não os implementava) e antes de eu ver T8; ao notar que T6 não correspondia à regra escrita, completei o pós-processamento em `analyze_synthetic_identifiability.py`.
    O resultado (`A_ONLY` 0/120 com H_A verdadeiro) foi uma surpresa de magnitude para mim (esperava poder reduzido, não nulo) e reforçou a recomendação; está declarado. O ramo L4 (largura) foi motivado por um princípio (o IC alargado tem de caber em 0,80 dex para separar A de B) mas formulado depois de eu conhecer as larguras de T2.

## 5. Inconsistências encontradas com os registros anteriores (nenhum arquivo existente foi editado)

| id | onde | o que diz | o que conflita | gravidade |
|---|---|---|---|---|
| I-1 | `PREREGISTRATION_STAGE2.md:3` | `Status: LOCKED` | o mesmo arquivo conserva texto de rascunho: título "(RASCUNHO)" (`:1`), "Não trava. `Status` permanece `DRAFT`" (`:549`), "vazio — este documento é um rascunho" (`:678`). Texto obsoleto num documento travado (herdado de `..._DRAFT.md`). | baixa (cosmética), mas um referee hostil pode usá-la |
| I-2 | `RESULTS_PRIMARY_STAGE2.md:45-50` | "**Dois** dos cinco gatilhos ... dispararam" e "Os outros **três** ... não dispararam" (lista quatro) | a própria tabela de `:192-200`, o debunker (`ADVERSARIAL_DEBUNKER_REPORT.md:122`) e `CLAIM_LEDGER.yaml:530` ("só o gatilho 3") dizem que **um** disparou (e quatro não) | baixa (contagem), mas está num resumo executivo |
| I-3 | `TEST_QUEUE.yaml:409-411` (e `CLAIM_LEDGER.yaml`, `known_gaps`: "holdout ... permanece disponivel para uma tentativa futura com calibracao de f_multi condicional a RUWE implementada") | holdout "disponível" para a tentativa condicional | `DISC-DEC-155`, `AGENTS.md` passo 9 e `03_REPLICATION_GATE/PROTOCOL.md`: o holdout só abre no Gate/lock próprio com autorização humana separada. "Disponível" pode ser lido como autorização. | média (ambiguidade de linha vermelha) |
| I-4 | `CLAIM_LEDGER.yaml` `DISC-CLAIM-008` (`counterevidence`, `exact_statement`) | o excesso do RUWE-alto no âncora é "evidência direta e não-ambígua de que o f_multi único não corrige adequadamente a heterogeneidade real" e que o `a0_fit` pooled "não pode ser interpretado como medição limpa" | **a simulação sintética deste rascunho reproduz o padrão do debunker (`D≈+0,28`) com `f` heterogêneo independente do resto, e nesse caso o `a0` escalar agrupado NÃO é enviesado** (cfg 0 vs 1: +0,213 vs +0,216 dex, p=0,88). A inferência "heterogeneidade ⇒ `a0` não confiável" só vale se a composição varia com `g_N` (cfg 10/11) ou se a *forma* da contaminação difere (cfg 3). O fechamento `CLOSED_INCONCLUSIVE` pode continuar justificado por outras razões (viés intrínseco, universalidade, GOF), mas a inferência específica não decorre **só** do achado RUWE. | média (lógica de uma inferência registrada). **Não é uma correção do registro** — é um ponto a levar ao autor/ao referee. |
| I-5 | `RESULTS_PRIMARY_STAGE2.md:246-249` | "Nenhum dado FÍSICO de nenhum sistema do holdout foi lido" | `run_stage2_primary_analysis.py:178-186,200`: `pd.read_parquet(QF_SAMPLE_PATH)` carrega o arquivo inteiro e `json.load` carrega o dicionário completo (com a lista de holdout), filtrando só depois; pela documentação (`PREREGISTRATION.md:94-100`) o parquet contém 43.147 linhas = 30.203 + 12.944. O grep de portão (`:137-141`) só verifica *subscrição* literal das chaves de holdout. Verdadeiro no sentido "não computado/usado", não "não carregado". | média (linha vermelha) |
| I-6 | `DIAGNOSTICO_AUTOMACAO_LINHAS_PESQUISA.md:61-66` | "SPARC-001 a SPARC-004 fecharam como INCONCLUSIVE devido à contaminação ... em estrelas com alto RUWE no Gaia DR3"; "fixar ... corte de RUWE < 1.4" | SPARC-001/003 fecharam por outras razões (sem observável discriminador; imagem do modelo estruturalmente fora do dado — `METHODOLOGY_EXTENSIONS.md` §1); SPARC-004 pelo RUWE. O corte de exclusão difere da calibração condicional, que é o caminho nomeado no `TEST_QUEUE.yaml`. (A nota datada de `DISC-DEC-155` em `:128-141` já trata da defasagem; estes dois detalhes não estão nela.) | baixa |
| I-7 | `DISC-DEC-155` ("seeds: use somente 20260955000-999") vs. `selfcal_pipeline.py:253-256,301,331,473` | seeds no bloco | o pipeline deriva sub-seeds por deslocamentos fixos que caem fora do bloco literal. Cumprível só no nível das seeds *passadas*. | baixa |
| I-8 | `RESULTS_SUMMARY_STAGE1.md:89-98` e `METHODOLOGY_ADDENDUM.md:239-260` vs. `PREREGISTRATION.md:8-18` | o viés residual `1,4-1,6×` do `a0` é atribuído ao modelo simplificado de wobble de fotocentro | o controle positivo de SPARC-004 já mostrava razão recuperado/previsto = **1,44** *sem* contaminação (`PREREGISTRATION.md:11-12`); **minha cfg 12 (`f=0`, nenhuma companheira, nenhum wobble) dá +0,212 ± 0,035 dex (fator 1,63)**, indistinguível da cfg 0/1 (p=0,96 / 0,86). O viés é, no mínimo, majoritariamente **independente de wobble**. Causa raiz não investigada por mim (hipótese a testar: o boost MOND da verdade é aplicado com `g_N` 3D, o ajuste usa `δ_AQUAL` com `g_N` projetado do bin). | **alta** (afeta a interpretação do Estágio 2 e a camada de viés) |
| I-9 | `PREREGISTRATION_STAGE2.md:35,380` (faixa `f_multi=0,25-0,47`) e `PROVENANCE_CHAE_EQS.md:213-225` (0,48/0,36) vs. `ADVERSARIAL_DEBUNKER_REPORT.md:192-195` ("69,2% ... RUWE<1,2") | a faixa é usada como referência para a amostra desta linha (gatilho 3) | não está documentado se as amostras de Chae que dão 0,36/0,48 aplicam corte de RUWE. Se aplicam, `f̂=0,10` (sem corte de RUWE) não é comparável. **Não verifiquei na fonte** (ver V3). | média (interpretação do gatilho 3) |
| I-10 | `PREREGISTRATION_STAGE2.md:250-294` (`bootstrap_a0_refit` sobre `final_raw` a `f̂` fixo) | IC95% de `a0` | não propaga a incerteza de `f̂` (o `f̂` não é recalibrado nas réplicas). Na simulação, o condicional com `f̂` fixo tem largura ≈ escalar (0,18-0,20 dex) e com recalibração ≈ 2× (0,33-0,40 dex, T2): o IC do Estágio 2 é, por construção, otimista na parte de calibração. Não afeta os vereditos já fechados (`CLOSED_INCONCLUSIVE`), mas afeta a leitura do IC publicado. | média |
| I-11 | procedimento desta sessão (não é do registro) | "Escreva APENAS em stage3_ruwe_conditional" | `pip download` criou transitoriamente um `.whl` em `COSMOLOGY_WIDE_BINARIES_MC_DEPROJECTION/`; removido imediatamente (D1). | baixa, mas é uma violação de escopo minha |

Consistentes (conferi): `DISC-DEC-155` ↔ `DISC-DEC-029/032` ↔ `TEST_QUEUE.yaml` ↔ `CLAIM_LEDGER.yaml` quanto a "19% / subgrupo RUWE-alto / CLOSED_INCONCLUSIVE / holdout intocado"; os números de `RESULTS_PRIMARY_STAGE2.md` batem com os do `TEST_QUEUE.yaml` e do `CLAIM_LEDGER.yaml`;
a diferença cru/corrigido do debunker (+0,237/+0,238) bate com a tabela; `5784/30203 = 0,19150`.

## 6. Riscos

R1. **"Inverse crime".** A verdade sintética e o estimador compartilham o mesmo modelo de contaminação; recuperação favorável aqui não implica recuperação no dado real. As mis-especificações testadas (3 botões) são as que eu escolhi.
R2. **Botão do artefato calibrado por semelhança qualitativa com os números publicados do debunker** (0,3 km/s). Isso torna o cenário "realista" circular em relação ao único dado real conhecido: **não** é evidência de que o artefato seja a causa, e o cenário foi escolhido depois de ver resultados de desenvolvimento (D2).
R3. **Poucas réplicas (10 por config).** Taxas como "10/10" têm limite inferior de Wilson de ~0,72; diferenças pequenas (p.ex. escalar vs. condicional em cfg 0/1) são não-significativas por construção, não prova de igualdade. As diferenças grandes (cfg 3, 10, 11) têm p < 10⁻⁵ (Welch).
R4. **Gates calibrados nas mesmas réplicas em que são avaliados** (30 nulas) — taxas nas configs 0, 1, 12 são in-sample. Limiares recalibrar com ≥ 200 réplicas.
R5. **A equivalência "agrupado = homogêneo" (S2) vale para rótulo de RUWE independente do resto.** No dado real, RUWE correlaciona com brilho, distância, separação etc. O fato real que decide é a composição por `g_N` (V2), que eu não pude ver.
R6. **A aritmética retrospectiva (S9) usa IC publicados de dois agentes diferentes** (N_boot=500 do debunker; 2000 do Estágio 2), Gaussianidade, sem covariâncias nem incerteza de `f̂`. O número `T=66` não deve ser citado como estatística de teste.
R7. **O viés intrínseco (+0,21 dex) é propriedade do estimador do pipeline travado** (confirmei por cfg 12 e pelo teste U1: nenhuma das minhas adições o causa), mas a causa raiz é desconhecida; se for um artefato do gerador sintético (e não do estimador sobre dado real), a interpretação de S4 mudaria. O controle positivo de SPARC-004 (razão 1,44, `PREREGISTRATION.md:11-12`) sugere que não é só meu gerador.
R8. **Protótipo com *monkeypatch***: funciona (U2), mas é frágil como implementação final (J12). O FPB usa interpolação linear em `f` (erro ≤0,004 dex em U3, num caso).
R9. **A recomendação "não travar" não é um veredito sobre física.** Ela é sobre a capacidade do desenho de produzir um veredito defensável. Um leitor apressado pode confundir; o rascunho repete a ressalva em §0, §5.4 e §10.
R10. **Extensão:** o rascunho tem ~600 linhas; conferi números contra `analysis_summary.log`/`synthetic_identifiability.log` e as citações `arquivo:linha` por leitura numerada, mas não por verificação automática de todas. Um referee deve checar por amostragem.
R11. **Nenhum referee hostil, nenhum segundo agente, nenhum spot-check independente** foram executados sobre este rascunho nem sobre o script (`DISC-DEC-155` os prevê depois de mim).
R12. **Declaração de "poder do gate" depende do jitter escolhido.** Para o artefato de 0,1 km/s, o poder é ≤ 4/10 e a regra completa depende então da largura, não dos gates (cfg 4: 5/10 `[GATE]`, 4/10 `[LARGURA]`).

## 7. Mapa para o critério de refutação da frente (DISC-DEC-155: (a)-(d))

- **(a) caminho para tocar holdout/dado real sem lock posterior.** Caminhos que considerei e como o rascunho os fecha: o portão de contagem E0 é o único acesso a dado real antes do resultado, só descoberta, só contagens, **depois** do LOCK (§4.3, §9.10);
  variantes de sensibilidade (K=3, bordas 1,3/1,5) só usam a descoberta (§6.10); o N do holdout aparece apenas como número em simulação (§8.4); `catalog.parquet` bruto exige autorização separada (§2.2, J16); linhas do holdout em memória via `read_parquet` completo é a **lacuna aberta** (§8.5, J14, I-5) — reconhecida, não resolvida.
- **(b) zona de veredito ambígua/reinterpretável pós-hoc.** A §5.2 define `V` como função determinística com ramos enumerados; §5.3 lista o que é `NAO_DISTINGUE` e o que é proibido mudar. **Fragilidades residuais:** parâmetros provisórios (`F`, `T*`, `zmax*`, `χ²*`, `N_A_MIN`) — a estrutura é pré-declarada, os números são julgamento aberto; L2 depende de definições operacionais de "borda" e "gatilho 9" que só ficam completas após J5/J6; o papel do debunker tem discricionariedade (não mecanizável).
- **(c) validação sintética incapaz de distinguir condicional de escalar.** Respondida em §0.3: distingue *heterogeneidade* (S1) e distingue o `a0` *quando a composição varia com `g_N`* (S3); **não** distingue o `a0` no caso-base (S2) — dito explicitamente como aplicação parcial do critério (c).
- **(d) inconsistência com `PREREGISTRATION_STAGE2.md`/`TEST_QUEUE.yaml`.** Inconsistências do **registro** em §5 acima. **Diferenças deliberadas** entre este rascunho e o Estágio 2 (cada uma justificada no texto; um referee pode vê-las como inconsistência):

| tema | Estágio 2 | este rascunho |
|---|---|---|
| modelo de contaminação | `f_multi` escalar | constante por partes em RUWE (K=2) |
| IC de `a0` | bootstrap com `f̂` fixo (`:250-294`) | FPB com recalibração (IC ≈ 2× mais largo) |
| ajuste de `a0` | 5 pontos, pesos iguais | 2×5 pontos, pesos = contagens |
| camada de viés | `a0`/IC ÷ 1,4-1,6, unilateral (`:333-353`) | `[lo/F, hi·F]`, `F=1,65`, bilateral |
| gatilho 4 | padrão agregado de `δ` | gate de universalidade entre estratos |
| gatilhos | só acionam o debunker | alguns (1, 2, 6, 7, 9, 10) também rebaixam |
| comparações | 2 hipóteses, 1 IC | 2 + 3 gates + 4 sensibilidades (declaradas) |
| linha vermelha | grep de subscrição | + varredura de dicionário, + extrato só-descoberta (J14) |

## 8. O que o autor humano precisa decidir primeiro (resumo operacional)

1. **Ler §0.4 do rascunho** e decidir entre: arquivar (J17c), reformular o estimador (J8) antes de qualquer lock, ou seguir só-diagnóstico (J17b).
2. **J1** (papel da descoberta) e **J8** (viés intrínseco) são pré-requisitos lógicos dos demais.
3. Nada neste diretório autoriza baixar dado, executar sobre dado real ou abrir o holdout. `Status` permanece `DRAFT`.
