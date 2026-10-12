# Histórico do desenvolvimento do GCM (Generalized Cell Mapping)

**Projeto:** IC FAPESP 2025/11057-0 — Ferramentas numéricas baseadas no método Cell-to-Cell para a obtenção de bacias de atração em sistemas dinâmicos não-lineares com um grau de liberdade
**Aluno:** Leonardo Franco Moreno · **Orientador:** Prof. Dr. Rodrigo Provasi Correia
**Repositório:** `leo-frm/PoliBoA-Leonardo-Moreno` (pastas `desenvolvimento_scm/` e `desenvolvimento_gcm/`)
**Documento gerado em:** 11/10/2026, como base para o relatório científico final

> **Atualização de 11/10/2026 (noite): decisão final do aluno.** SCM e GCM devem ficar o mais parecidos possível com o PoliBoA original.
> - **O integrador voltou a ser o do original** (padrão do DifferentialEquations, sem método e sem tolerância). O Vern9 1e-8 descrito em várias seções abaixo foi **revertido**.
> - **Ficam:** o agrupamento por trajetória contínua (Hsu §8.2), apresentado como parte do método SCM, e o `reset_dt` do GCM.
> - **A M5 foi enxugada** e os diagnostics do GCM foram apagados.
> - **Onde ver:**
>   - novos números: §2.5, em "Decisão final";
>   - mudanças: §5, §6, §7 e §8;
>   - estudo de precisão com Vern9: Anexo A.
>
> As §§1–4 foram mantidas como registro histórico da rodada com Vern9.

---

## Como ler este documento

Todo número vem acompanhado da sua fonte entre colchetes. As fontes usadas são:

| Marcador | O que é | Onde está |
|---|---|---|
| `[git:<hash>]` | Commit do repositório | `git log` / `git show <hash>` |
| `[ctx:CLAUDE]` | Resumo do estado do projeto preparado antes desta sessão | `poliboa_contexto_claude_code/CLAUDE.md` |
| `[ctx:HIST]` | Resumo das conversas jan–out/2026 | `poliboa_contexto_claude_code/docs/HISTORICO.md` |
| `[ctx:REL]` | Texto do relatório parcial | `poliboa_contexto_claude_code/docs/relatorio_parcial.md` |
| `[ctx:PROTO]` | Log do protótipo Python (experimento do orientador) | `poliboa_contexto_claude_code/docs/prototipo_python/exp_hd_1000.log` |
| `[term:10/10]` | Saída de terminal do aluno, lida durante a sessão de 10/10/2026 | não foi salva em arquivo (registro da conversa) |
| `[term:11/10]` | Idem, sessão de 11/10/2026 | idem |
| `[diag:ref]` | Execução de `desenvolvimento_gcm/diagnostic_ref_gcm.jl` (10/10) | saída de terminal (registro da conversa); figuras `diag_ref_*.png` |
| `[diag:tol]` | Execução de `desenvolvimento_gcm/diagnostic_tol_gcm.jl` (10/10) | saída de terminal (registro da conversa); figura `diag_tol_scm_hd_1000.png` |
| `[sp:s1]` … `[sp:s5]` | Scripts de investigação executados por Claude (10/10) | pasta temporária `SCRATCH` (ver abaixo); só o s5 tem log em arquivo (`s5.log`) |
| `[log:m4]`, `[log:m5]` | Execução oficial de `tests_milestone4_gcm.jl` e `tests_milestone5_gcm.jl` (11/10) | `SCRATCH/m4.log`, `SCRATCH/m5.log` |
| `[img:aluno-1]`, `[img:aluno-2]` | Imagens enviadas pelo aluno no chat (10/10) | `C:\Users\Leonardo\AppData\Local\Temp\claude\…\images\1.png` e `2.png` |
| `[conv]` | Afirmação que só existe na conversa, sem arquivo | — |

`SCRATCH` = `C:\Users\Leonardo\AppData\Local\Temp\claude\C--Users-Leonardo-Desktop-PoliBoA---Leonardo-Moreno\c22e2f95-e3a6-4f00-952f-d028ef7536c7\scratchpad\`

> ⚠ **A pasta `SCRATCH` é temporária** e pode ser apagada pelo sistema. Se os scripts e logs forem citados no relatório, copie-os para um lugar permanente (fora do GitHub, se preferir).

Convenções de interpretação:

- **⚠ Interpretação:** conclusão minha, apoiada nos dados mas não medida diretamente.
- **⚠ Incerteza:** dado não confirmado ou com fonte frágil.
- **(derivado):** número calculado a partir de outros números citados, com a conta indicada.

---

# Parte 1 — Visão geral e linha do tempo

## 1. Visão geral

### 1.1 Objetivo

Estender o módulo de Cell Mapping do PoliBoA do **Simple Cell Mapping (SCM)** para o **Generalized Cell Mapping (GCM)** de Hsu (1987, caps. 10–12), mantendo a mesma estrutura de quatro milestones do SCM:

- **M1:** construção do mapa;
- **M2:** detecção de atratores e bacias;
- **M4:** paralelização por partição;
- **M5:** validação de tempo, memória e precisão.

A M3 (lazy evaluation) foi descartada em 01/04 `[ctx:HIST]`.

A diferença central entre os dois métodos:

- **SCM:** cada célula tem **uma** imagem, a célula onde cai o centro dela depois de um período. O mapa é um grafo com N arestas.
- **GCM:** cada célula tem **várias** imagens, com probabilidades. Integram-se várias amostras dentro da célula, e a fração de amostras que cai em cada célula-destino é a probabilidade de transição. O mapa vira uma **cadeia de Markov**:
  - atratores são **grupos persistentes** (conjuntos fechados de células que se comunicam);
  - cada célula tem uma **probabilidade de absorção** em cada atrator;
  - células de fronteira podem ter mais de um destino possível (**domicílio múltiplo**).

### 1.2 Ponto de partida

- **SCM concluído e entregue no relatório parcial** `[ctx:REL]`:
  - concordância de 86,8% com o PoliBoA original (HD 200×222);
  - memória de ≈ 0,66 MB/kCell, constante com o tamanho do grid;
  - speedup de 4,40× com 6 threads.
- **Problemas conhecidos do SCM** `[ctx:REL]`:
  - um atrator duplicado no Helmholtz-Duffing (HD), atribuído no relatório a "variações numéricas entre arquiteturas de hardware";
  - no Spring Column (SC), 9 atratores com 10 períodos de transiente, corrigidos com 100 períodos.
- **Primeiro commit do GCM** `[git:9de06a1]`: M1 e M2 em `desenvolvimento_gcm/`, com a M2 ainda como "SCM dominante" (ver §2.2).

### 1.3 Estado final (11/10/2026)

- **GCM completo nas quatro milestones**, seguindo a teoria de Hsu como cadeia de Markov, com amostragem de m subcélulas sorteadas sem reposição (proposta do orientador).
- **Testes oficiais com o código final** `[log:m4]`, `[log:m5]`, `[term:11/10]`:
  - HD com **4 atratores**, SC com **2**;
  - reprodutibilidade exata entre 1 e 4 threads (0 células diferentes);
  - memória linear em N (≈ 11,7 MB/kCell de alocação total no HD);
  - speedup de 3,39× (HD) e 3,26× (SC) com 4 threads.
- **Estudo de refinamento proposto pelo orientador (Teste 6), HD:**
  - um GCM grosso 200×100 com m = 5 (10% do custo do fino 1000×1000) erra a fração de área das bacias em **≤ 0,38 p.p.**;
  - erro médio por célula de 0,095;
  - esses números coincidem com o protótipo Python `[log:m5]`, `[ctx:PROTO]`.
- **Problema grave encontrado e resolvido nesta fase: atratores duplicados em grids finos.** A causa é uma propriedade conhecida do SCM (Hsu §8.2), e afetava também o SCM entregue. A correção vale para SCM e GCM (§2.5).
- **Integração numérica:** passou do padrão do DifferentialEquations (reltol 1e-3) para **Vern9 com reltol = abstol = 1e-8** no SCM e no GCM. O PoliBoA original foi mantido intocado como **grupo controle** (decisão do aluno).

---

## 2. Linha do tempo das milestones

### 2.0 Cronologia resumida

| Quando | Versão | O que aconteceu | Fonte |
|---|---|---|---|
| jan–ago/2026 | SCM | M1, M2, M4 e M5 do SCM; relatório parcial | `[ctx:REL]`, `[git]` |
| 15/08/2026 | SCM | `compare_basins` corrigido (guloso bijetivo, −1 ↔ −1 fixo); `compute_single_cell_mapping` movido para `utils.jl` | `[git:15f367b]` |
| 15/08/2026 | **GCM v0** | Criação de `desenvolvimento_gcm/`: estruturas CSR, amostragem em retículo S^d, M1, e M2 = SCM com "imagem dominante" | `[git:9de06a1]`, `[git:fdc51bc]` |
| set/2026 (⚠ data incerta) | **GCM v1** | M2 reescrita pela teoria de Hsu (cadeia de Markov); M4 e M5 criadas; retículo S×S com S = 3 | `[ctx:CLAUDE]`, `[ctx:HIST]`; branch remota `gcm-milestones-2-4-5` (commit 5e48742, não inspecionado localmente) |
| 10/10/2026 | — | Reunião com o orientador: sortear m de K subcélulas; comparar GCM grosso × SCM refinado | `[ctx:CLAUDE]` |
| 10/10/2026 | **GCM v2** | Sorteio sem reposição, `reset_dt` no integrador, Testes 0 e 6 na M5 (escrita no claude.ai, sem rodar Julia) | `[ctx:CLAUDE]` |
| 10/10/2026 | v2 rodada | Primeiras execuções em Julia: M1, M4 e M5 passam, mas a referência do Teste 6 tem 5 atratores | `[term:10/10]` |
| 10/10/2026 | **v2.1** | Investigação dos atratores duplicados: integrador → Vern9 1e-8; agrupamento por trajetória contínua (Hsu §8.2) no SCM e no GCM | `[diag:*]`, `[sp:s1–s5]` |
| 11/10/2026 | v2.1 | Testes oficiais M1, M4 e M5 com o código final | `[term:11/10]`, `[log:m4]`, `[log:m5]` |

> ⚠ **Incerteza de datas:** o `HISTORICO.md` situa a criação do GCM entre 27/09 e 10/10 e a correção do `compare_basins` em 09/08. O git registra as duas em 15/08/2026 (`9de06a1`, `15f367b`). Vale conferir qual data usar no relatório.

---

### 2.1 M1 — Construção do mapa (`milestone1_gcm.jl`, `structures_gcm.jl`, `utils_gcm.jl`)

**Contexto.** No SCM, `compute_single_cell_mapping` integra **o centro** de cada célula por (transiente + 1) períodos e registra a célula de destino. No GCM, é preciso integrar **várias amostras** por célula e guardar todos os destinos com suas frequências.

**Motivação / hipótese.** Hsu (§11.1) representa o GCM por três arranjos:

- **I(z):** número de imagens da célula z;
- **C(z,i):** a i-ésima imagem;
- **P(z,i):** a probabilidade de transição para essa imagem.

Guardá-los em formato **CSR** (listas concatenadas com deslocamentos) mantém a memória proporcional a N·Ī, em que Ī é o número médio de imagens por célula. Isso preserva o argumento de memória O(N) do relatório parcial.

**Alterações, por versão:**

- **v0** `[git:9de06a1]`:
  - `GcmProblem(base, samples_per_side)`, com S amostras por dimensão, formando um **retículo regular** de S^d pontos na célula;
  - `GeneralizedCellMap` em CSR (`offsets`, `images`, `counts` como `Int32`, `computed`);
  - `store_cell_sample!` (retículo, Hsu §11.12);
  - `map_point`, que generaliza `compute_single_cell_mapping` para um ponto qualquer;
  - `build_gcmap` e `compute_cell_transition_probabilities`;
  - o destino −1 (célula sumidouro) passa a ser **uma imagem legítima**, e não um descarte.
- **v1** `[ctx:CLAUDE]`: mesmo esquema, com retículo S×S e S = 3 (9 amostras por célula).
- **v2** `[ctx:CLAUDE]`, proposta do orientador em 10/10:
  - `GcmProblem(base, subgrid, m; seed = 1, sampling = :random | :stratified)`. A célula é dividida em K = prod(subgrid) subcélulas (ex.: 10×10 = 100), e **sorteiam-se m delas sem reposição** (ex.: m = 10).
  - `select_samples(gp, cell_id)`: o gerador aleatório é `Xoshiro(seed·10⁹ + cell_id)`. Assim o sorteio é **reprodutível e independente da ordem de processamento e do número de threads**. Com `:stratified`, faz-se um sorteio em cada uma de m faixas contíguas de subcélulas.
  - `store_cell_sample!` grava o **centro da subcélula**. Quando o grid fino de referência é o grid grosso × subgrid, cada amostra coincide com o centro de uma célula fina; o Teste 6 explora isso.
  - `map_point` passa a reiniciar o integrador por completo a cada ponto: `reinit!(…; reset_dt = true)`. Motivação: o problema de 1 thread ≠ 4 threads da v1 (ver M4).
- **v2.1** (esta sessão, 10/10): `make_gcm_integrator` usa **`Vern9()` com `reltol = abstol = 1e-8`**, no lugar do padrão do DifferentialEquations (reltol 1e-3, abstol 1e-6). Motivação e evidência na §2.5.

**Testes / resultados** (`tests_milestone1_gcm.jl`: HD 200×222 e SC 200×200, subgrid 10×10, m = 10):

| Versão | HD: atratores e bacias (células) | SC: atratores e bacias (células) | Fonte |
|---|---|---|---|
| v1 (S = 3) | 4: 25354 / 5228 / 9340 / 4478, nenhuma divergente | 2: 14134 / 25638, 228 divergentes | `[ctx:CLAUDE]` |
| v2, tolerância padrão | 4: 4947 / 5425 / 9332 / 24696 (soma 44400, 0 divergentes) | 2: 15495 / 24212; 293 divergentes (derivado: 40000 − 15495 − 24212) | `[term:10/10]` |
| v2.1 (Vern9 + agrupamento) | 4: 4955 / 24530 / 9346 / 5569 | 2: 25719 / 13815; 466 divergentes (derivado: 40000 − 25719 − 13815) | `[term:11/10]` |

Sobre a última linha:
- O número da bacia verde (5569) ficou ilegível no terminal. Ele foi derivado como 44400 − 4955 − 24530 − 9346. A execução equivalente em `[sp:s5]` mostra que não há células divergentes no HD (frações 0,1116 / 0,5525 / 0,2105 / 0,1254).
- ⚠ **Interpretação:** as 466 células divergentes do SC (≈ 1,2%) vêm da "faixa morta" da região [−3,14; 3,14] (ver §4.6), e não de divergência real.

**Decisão.** A M1 v2.1 é a versão final. O sorteio por semente de célula garantiu a reprodutibilidade (ver M4).

---

### 2.2 M2 — Detecção de atratores e bacias (`milestone2_gcm.jl`)

**Contexto.** No SCM, `find_attractors_from_scm` segue o ponteiro único de cada célula até:

- encontrar uma célula já classificada (a trajetória herda o rótulo dela);
- sair do grid (vira −1);
- ou repetir células `attractor_tolerance` vezes (atrator novo).

**v0, abandonada** `[git:9de06a1]`. A M2 copiava integralmente o algoritmo do SCM, trocando o ponteiro `scm.target[current]` pela **imagem de maior probabilidade** (`dominant_image`).

**Por que foi abandonada** `[ctx:HIST]`:
- Usar só a imagem dominante transforma o GCM no **"SCM dominante"** de Hsu (§11.3): um mapa determinístico selecionado a partir do GCM, que joga fora as probabilidades.
- Decisão registrada: *"o GCM correto não é determinístico (não é SCM dominante)"*.
- Além disso, essa mudança de mecanismo foi feita **sem aviso** pelo assistente na época, o que gerou a regra "não mudar algoritmo sem avisar antes".

**v1/v2, teoria de Hsu completa** (caps. 10–11). Fluxo de `find_attractors_from_gcm` → `analyze_gcm` → `basin_result_from_analysis`:

1. **`build_preimage_array`** (§11.2): arranjo de pré-imagens em CSR, ou seja, quem mapeia para cada célula. O sumidouro é o nó N+1 (convenção de Hsu §11.1).
2. **`analyze_dominant_scm`** (§11.3 e §8.2): análise do SCM dominante, que fornece as **soluções periódicas** usadas como sementes. Pelo Teorema 11.3.1, todo grupo persistente contém ao menos uma.
3. **`identify_persistent_groups`** (§11.4–11.5). Para cada semente:
   - expande-se o conjunto alcançável seguindo **todas** as imagens;
   - se ele toca uma célula já classificada, a semente é transiente;
   - se é fechado e todas as suas células voltam à semente (busca para trás), é um **grupo persistente**.
   - É equivalente a encontrar as componentes fortemente conexas fechadas do grafo. O sumidouro é o grupo 1.
4. **`group_period`** (§11.6): período pelo máximo divisor comum das diferenças de nível em uma busca em largura.
5. **`classify_markov_group_type`**: regular ou quase-periódico, com o mesmo corte de 20 células do SCM.
6. **`determine_domiciles`** (§11.10): busca para trás a partir de cada grupo. Dm(z) é o número de grupos alcançáveis a partir de z.
7. **`assign_basin_probabilities`** (§11.9 e §11.11.2):
   - para células persistentes ou de domicílio único, a absorção é trivial (probabilidade 1 no único grupo possível);
   - para as de domicílio múltiplo, resolve-se o sistema esparso **(I − Q)X = B**.
8. **`basin_result_from_analysis`**: o rótulo de cada célula é o argmax da probabilidade de absorção; o grupo 1 vira −1. A saída é um `BasinResult` comum, o que permite reaproveitar `plotbasins`, `report_attractors` e `compare_basins` do SCM.

**v2.1, revisão e correção (10/10):**

- **Revisão da lógica contra Hsu** (sem alterar código): grupos persistentes, período, domicílios e absorção **estão corretos**. Pontos encontrados:
  - empate exato no argmax manda a célula para o sumidouro (gravidade baixa, não alterado);
  - o rótulo da referência do Teste 6 estava enganoso (corrigido, §2.4);
  - comparações GCM × SCM misturam diferenças de integrador (resolvido com Vern9 nos dois).
- **Nova etapa `merge_clustered_groups`**, logo após `identify_persistent_groups`. Para cada grupo persistente, integra-se a **trajetória contínua** (sem reiniciar no centro da célula); se ela termina numa célula de outro grupo, os dois são o mesmo atrator e são unidos. Motivação na §2.5. Com m = 1 (equivalente ao SCM), isso é necessário. Com m = 10, as amostras normalmente já conectam as células, e a etapa não muda nada `[sp:s2]`.

**Resultados, análise do GCM** (Teste 2 da M5, HD 200×222 e SC 200×200, m = 10):

| | HD v2 (tol. padrão) | HD v2.1 | SC v2 (tol. padrão) | SC v2.1 |
|---|---|---|---|---|
| Grupos persistentes (com sumidouro) | 5 | 5 | 3 | 3 |
| Células persistentes | 7 (0,02%) | 4 (0,01%) | 3 (0,01%) | 3 (0,01%) |
| Domicílio único | 63,83% | 64,61% | 39,24% | 72,77% |
| Domicílio múltiplo (fronteira) | 36,17% | 35,39% | 60,76% | 27,23% |
| Fonte | `[term:10/10]` | `[log:m5]` | `[term:10/10]` | `[log:m5]` |

⚠ **Interpretação:** a queda do domicílio múltiplo no SC (60,76% → 27,23%) mostra que boa parte da "fronteira" do SC era **ruído de integração**. Com reltol 1e-3 e 100 períodos de transiente, amostras vizinhas caíam em destinos diferentes por erro numérico, e não pela dinâmica.

**Decisão.** A M2 v2.1 é a versão final. A etapa de agrupamento é a única adição em relação à teoria dos caps. 10–11, e ela também é fundamentada em Hsu (§8.2).

---

### 2.3 M4 — Paralelização por partição (`milestone4_gcm.jl`)

**Contexto.** No SCM, cada thread processa uma faixa contígua de células com um integrador próprio, e os mapas parciais são unidos em série (`join_scmaps!`). O GCM segue a mesma ideia.

**Alterações (v1/v2):**
- `process_cell_range_gcm`: cada thread processa uma faixa e constrói um CSR parcial.
- `join_gcmaps!`: concatena os CSR parciais **em ordem**, deslocando os offsets. É diferente do `join_scmaps!`, que preenchia posições vazias em qualquer ordem; aqui a ordem de união importa, e o `foldl` a preserva.
- `build_gcmap_parallel`: mesma interface do `build_scmap_parallel`.

**Problema da v1** `[ctx:CLAUDE]`. 1 thread ≠ 4 threads: mais células no grupo do atrator verde.
- Hipótese: o integrador guardava o passo adaptativo (dt) da integração anterior, então o resultado dependia de **qual célula foi integrada antes**, ou seja, da partição.
- Correção (v2): `reinit!(…; reset_dt = true)` em `map_point`, somada ao sorteio por semente de célula.

**Testes / resultados:**

| Teste | Resultado | Fonte |
|---|---|---|
| Teste 0 (M5), v2 tol. padrão, HD | 0 de 44400 células com mapa diferente; `compare_basins` = 1,0000 | `[term:10/10]` |
| Teste 0 (M5), v2 tol. padrão, SC | 0 de 40000; 1,0000 | `[term:10/10]` |
| Imagens M1 (1 thread) × M4 (4 threads), v2 | idênticas byte a byte (HD e SC) | comparação com `cmp` na sessão de 10/10 `[conv]` |
| Teste 0 (M5), v2.1, HD e SC | 0 de 44400 e 0 de 40000; 1,0000 | `[log:m5]` |
| Imagens M1 × M4, v2.1 | idênticas byte a byte (HD e SC) | `cmp` em 11/10; `gcm_basin_*_testm{1,4}.png` |
| Tempo da execução da M4, v2.1 | 7 min 14 s | `[log:m4]` (saída do `time`) |

**Escalabilidade por threads** (Teste 5 da M5, HD 200×222 e SC 200×200, m = 10):

| Threads | HD v2 (padrão) | HD v2.1 | SC v2 (padrão) | SC v2.1 |
|---|---|---|---|---|
| 1 | 15,36 s (1,00×) | 42,53 s (1,00×) | 86,65 s (1,00×) | 247,35 s (1,00×) |
| 2 | 7,72 s (1,99×) | 20,83 s (2,04×) | 43,01 s (2,01×) | 124,33 s (1,99×) |
| 4 | 4,73 s (3,25×) | 12,55 s (3,39×) | 26,64 s (3,25×) | 75,95 s (3,26×) |
| Fonte | `[term:10/10]` | `[log:m5]` | `[term:10/10]` | `[log:m5]` |

⚠ Esses números não se comparam diretamente com o 4,40× do SCM no relatório parcial, que foi medido com 6 threads `[ctx:REL]`.

**Decisão.** A hipótese do dt foi **confirmada**: com `reset_dt`, o resultado não depende do número de threads. A M4 é a versão final.

---

### 2.4 M5 — Validação (`milestone5_gcm.jl`, `tests_milestone5_gcm.jl`)

> ⚠ Os resultados desta seção são das rodadas de 10/10 (tolerância padrão) e de 11/10 de manhã (Vern9), com a M5 antiga (Teste 0 e Teste 6 também estratificado). A M5 foi enxugada depois, e os resultados finais estão na §2.5, em "Decisão final".

**Funções** (sufixo `_gcm` para as que diferem do SCM):

- `benchmark_method_gcm`: conta as integrações como N·m.
- `memory_scaling_gcm`: reporta também Ī e a memória do próprio mapa.
- `thread_scaling_gcm`.
- `basin_id_mapping` e `compare_basins_by_domicile`: separam a concordância entre células de domicílio único e múltiplo.
- `reproducibility_gcm` (Teste 0).
- `refinement_study_gcm` (Teste 6).
- Reaproveitadas do SCM sem alteração: `BenchmarkResult`, `compare_basins`, `compare_benchmarks`.

Pela regra do projeto, a M5 é o **único** lugar com comparações entre métodos (GCM × SCM × original).

**Alterações nesta sessão:**
- `tests_milestone5_gcm.jl`: **só a ordem de execução mudou**. Agora roda Teste 0 (HD e SC) → Testes 1–5 (HD, depois SC) → Teste 6 por último, que é o mais caro. O cabeçalho foi atualizado.
- `milestone5_gcm.jl`: só o texto impresso da referência do Teste 6, que agora diz *"SCM refinado …, calculado como GCM com m = 1"*.

#### Resultados dos Testes 1–5 (HD 200×222, SC 200×200, subgrid 10×10, m = 10)

**Teste 1, benchmark do GCM** (3 repetições):

| | HD v2 (padrão) | HD v2.1 | SC v2 (padrão) | SC v2.1 |
|---|---|---|---|---|
| Tempo médio | 14,88 ± 0,09 s | 41,90 ± 0,05 s | 86,42 ± 0,55 s | 249,93 ± 0,55 s |
| Bytes alocados | 546,28 MB | 546,04 MB | 480,73 MB | 472,25 MB |
| Integrações | 444000 | 444000 | 400000 | 400000 |
| Fonte | `[term:10/10]` | `[log:m5]` | `[term:10/10]` | `[log:m5]` |

**Teste 2, concordância de bacias:**

| | HD v2 | HD v2.1 | SC v2 | SC v2.1 |
|---|---|---|---|---|
| GCM × Original (grid inteiro) | 84,73% | 83,43% | 26,91% | 22,85% |
| GCM × Original (domicílio único) | 99,72% | 98,97% | 31,73% | 23,99% |
| GCM × SCM (grid inteiro) | 85,08% | 86,76% | 74,92% | **95,71%** |
| GCM × SCM (domicílio único) | 99,75% | 99,97% | 96,62% | 99,47% |
| GCM × SCM (domicílio múltiplo) | 59,20% | 62,63% | 60,91% | 85,67% |
| Fonte | `[term:10/10]` | `[log:m5]` | `[term:10/10]` | `[log:m5]` |

Como ler esta tabela:
- ⚠ **Interpretação:** no domicílio único, a concordância com o SCM fica entre 99,5% e 100%. As discordâncias se concentram na fronteira (domicílio múltiplo), como Hsu prevê (§10.6).
- ⚠ **Interpretação:** a concordância com o **original** fica limitada a ~83% no HD por causa da imprecisão do **próprio original** (reltol 1e-3). Quando integrado com Vern9 (numa cópia), o original muda 16,6% das suas células, e o SCM concorda 98,3% com ele (§2.5, `[sp:s4]`).
- No SC, o original falha por outro motivo (§4.6).

**Teste 3, razões de custo, v2.1** `[log:m5]`:

| | Tempo | Bytes alocados | Integrações |
|---|---|---|---|
| HD: GCM / Original | 42,94 s / 1,35 s = 31,82× | 22,90× | 10,00× |
| HD: GCM / SCM | 42,94 s / 4,27 s = 10,06× | 17,64× | 10,00× |
| SC: GCM / Original | 250,59 s / 8,39 s = 29,85× | 20,83× | 10,00× |
| SC: GCM / SCM | 250,59 s / 24,55 s = 10,21× | 16,89× | 10,00× |

Para comparação, com a tolerância padrão (v2) as razões eram `[term:10/10]`:
- GCM / Original: 11,21× (HD) e 9,90× (SC);
- GCM / SCM: 10,56× (HD) e 10,30× (SC).

A razão GCM/SCM ≈ 10× bate com a razão de integrações, m = 10. A razão contra o original subiu porque **o original continua com reltol 1e-3** (grupo controle), enquanto SCM e GCM passaram a usar Vern9 1e-8.

**Teste 4, escalabilidade de memória:**

| Grid | MB/kCell HD v2 | MB/kCell HD v2.1 | Ī HD v2 → v2.1 | Mapa (MB) HD v2 → v2.1 | Fonte |
|---|---|---|---|---|---|
| 100×111 | 11,806 | 11,813 | 5,225 → 5,119 | 0,486 → 0,477 | `[term:10/10]`, `[log:m5]` |
| 200×222 | 11,734 | 11,728 | 5,300 → 5,152 | 1,970 → 1,920 | idem |
| 400×444 | 11,674 | 11,663 | 5,390 → 5,170 | 8,001 → 7,705 | idem |

| Grid | MB/kCell SC v2 | MB/kCell SC v2.1 | Ī SC v2 → v2.1 | Mapa (MB) SC v2 → v2.1 | Fonte |
|---|---|---|---|---|---|
| 100×100 | 11,453 | 11,316 | 3,547 → 2,618 | 0,310 → 0,239 | idem |
| 200×200 | 11,461 | 11,259 | 4,603 → 2,976 | 1,562 → 1,065 | idem |
| 400×400 | 11,536 | 11,237 | 5,226 → 2,998 | 7,009 → 4,289 | idem |

Leitura:
- MB/kCell constante significa que a alocação total é **O(N)**. Esse valor mede tudo o que é alocado na execução (dominado pela integração), não o tamanho do mapa.
- O mapa em si ocupa poucos MB: 7,7 MB para 177600 células.
- ⚠ **Interpretação:** com a tolerância padrão, o **Ī do SC crescia com N** (3,5 → 5,2), o que ameaçava a linearidade da memória do mapa. Com Vern9, Ī fica estável (2,6 → 3,0). O crescimento era artefato numérico.

**Teste 5:** ver a tabela de threads na §2.3.

#### Teste 6, estudo de refinamento (proposta do orientador, 10/10)

**Pergunta.** Um GCM grosso, com m amostras por célula, chega a resultados comparáveis aos de um SCM refinado, gastando menos?

**Configuração:**
- Referência: SCM 1000×1000, calculado como GCM com subgrid [1,1] e m = 1, para usar o mesmo integrador.
- GCMs grossos: 100×50 (subgrid 10×20) e 200×100 (subgrid 5×10), alinhados ao fino. Cada amostra do grosso é o centro de uma célula fina, então as diferenças vêm do **método**, não da integração.

**Métricas:**
- **Atratores:** quantos a referência tem e quantos o grosso encontra.
- **Erro fração:** maior diferença na área de uma bacia, em pontos percentuais. A área estimada vem da **probabilidade de absorção**, não do argmax.
- **Erro/célula:** distância de variação total média entre a absorção estimada e a fração real das subcélulas.
- **Ponto a ponto:** cada célula fina recebe o rótulo da sua célula grossa.

> **Atenção (fácil de confundir no relatório):**
> - O "ponto a ponto" do Teste 6 compara **resoluções diferentes**. Ele tem um teto: uma célula grossa 200×100 cobre 50 células finas e recebe um único rótulo, então nas regiões fractais ela necessariamente discorda de parte delas. Mesmo integrando **todas** as subcélulas (m = K), o ponto a ponto fica em 83,22% (200×100) e 73,37% (100×50) `[log:m5]`.
> - Isso não é a mesma coisa que a concordância SCM × original (~98%), que compara **o mesmo grid**.

**Resultados HD, v2.1** `[log:m5]`. Referência com **4 atratores**, frações 0,2090 / 0,1233 / 0,5458 / 0,1219.

| Grid | m | Integrações (custo) | Sorteio | Atratores | Erro fração | Erro/célula | Ponto a ponto |
|---|---|---|---|---|---|---|---|
| 100×50 | 5 | 25000 (2,5%) | random | 4 de 4 | 1,48 p.p. | 0,1226 | 76,77% |
| 100×50 | 10 | 50000 (5%) | random | 3 de 4 | 12,33 p.p. | 0,1595 | 72,15% |
| 100×50 | 20 | 100000 (10%) | random | 3 de 4 | 12,33 p.p. | 0,1476 | 72,62% |
| 100×50 | 200 | 1000000 (100%) | random | 3 de 4 | 12,33 p.p. | 0,1234 | 73,37% |
| 200×100 | 5 | 100000 (10%) | random | 4 de 4 | **0,38 p.p.** | 0,0948 | 81,01% |
| 200×100 | 10 | 200000 (20%) | random | 4 de 4 | 0,35 p.p. | 0,0631 | 82,22% |
| 200×100 | 20 | 400000 (40%) | random | 4 de 4 | 0,32 p.p. | 0,0402 | 82,77% |
| 200×100 | 50 | 1000000 (100%) | random | 4 de 4 | 0,21 p.p. | 0,0094 | 83,22% |
| 100×50 | 5 | 25000 | stratified | 3 de 4 | 12,33 p.p. | 0,1879 | 71,57% |
| 100×50 | 10 | 50000 | stratified | 3 de 4 | 12,33 p.p. | 0,1569 | 72,62% |
| 100×50 | 20 | 100000 | stratified | 3 de 4 | 12,19 p.p. | 0,1630 | 73,19% |
| 200×100 | 5 | 100000 | stratified | 4 de 4 | 0,36 p.p. | **0,0827** | 81,72% |
| 200×100 | 10 | 200000 | stratified | 4 de 4 | 0,28 p.p. | 0,0537 | 82,52% |
| 200×100 | 20 | 400000 | stratified | 4 de 4 | 0,29 p.p. | 0,0329 | 82,97% |

**Comparação com o protótipo Python** `[ctx:PROTO]`. Ele usa RK4 de passo fixo com 200 passos por período, bem mais preciso que reltol 1e-3; isso aparece no código em `gcm_proto.py`, `rk4_period(…, nsteps=200)`.

| | Protótipo | Julia v2.1 |
|---|---|---|
| Frações da referência | 0,209 / 0,1233 / 0,5458 / 0,1219 | 0,2090 / 0,1233 / 0,5458 / 0,1219 |
| 200×100, m = 5: erro/célula; ponto a ponto | 0,0954; 81,18% | 0,0948; 81,01% |
| 200×100, m = 10 | 0,0640; 82,19% | 0,0631; 82,22% |
| 200×100, m = 20 | 0,0401; 82,80% | 0,0402; 82,77% |
| 200×100, m = 50 | 0,0094; 83,21% | 0,0094; 83,22% |
| 100×50 | perde 1 atrator para m = 5, 10, 20 e 200 | perde 1 atrator para m = 10, 20 e 200; com m = 5 random encontrou os 4 |

Leitura dos resultados:
- ⚠ **Interpretação (100×50):** o resultado depende da semente; o protótipo já registrava que, no 100×50 com sorteio aleatório, o número de atratores varia com a semente `[ctx:CLAUDE]`. O atrator perdido tem a bacia de 0,1233 (o "amarelo", em (−0,74; −0,49)), e o erro de 12,33 p.p. é exatamente essa área.
- **Estratificada × aleatória** (derivado da tabela): o erro por célula cai 12,8% com m = 5 (0,0827 contra 0,0948), 14,9% com m = 10 e 18,2% com m = 20. O protótipo indicava ≈ 13% `[ctx:CLAUDE]`.

**Resultados SC, v2.1** `[log:m5]`. Referência com **2 atratores**, frações −1: 0,0260 / 1: 0,3510 / 2: 0,6230.

| Grid | m | Custo | Atratores | Erro fração | Erro/célula | Ponto a ponto |
|---|---|---|---|---|---|---|
| 100×50 | 5 | 2,5% | 2 de 2 | 0,11 p.p. | 0,1019 | 80,68% |
| 100×50 | 10 | 5% | 2 de 2 | 0,08 p.p. | 0,0716 | 81,98% |
| 100×50 | 20 | 10% | 2 de 2 | 0,32 p.p. | 0,0503 | 82,55% |
| 100×50 | 200 | 100% | 2 de 2 | 0,00 p.p. | 0,0000 | 83,21% |
| 200×100 | 5 | 10% | 2 de 2 | 0,04 p.p. | 0,0670 | 87,81% |
| 200×100 | 10 | 20% | 2 de 2 | 0,02 p.p. | 0,0453 | 88,44% |
| 200×100 | 20 | 40% | 2 de 2 | 0,03 p.p. | 0,0278 | 88,92% |
| 200×100 | 50 | 100% | 2 de 2 | 0,00 p.p. | 0,0000 | 89,17% |

**Primeira execução do Teste 6 (v2, tolerância padrão), descartada** `[term:10/10]`:
- A referência HD tinha **5 atratores** (frações 0,5528 / 0,1237 / 0,0244 / 0,0918 / 0,2073).
- O erro fração de 2,44 p.p. do 200×100 era exatamente a bacia espúria de 0,0244.
- Esse resultado levou à investigação da §2.5.

**Decisão.** O Teste 6 v2.1 é o resultado a levar ao orientador.
- No HD, um GCM 200×100 com m = 5 (10% do custo do fino) reproduz as frações de área com erro ≤ 0,38 p.p.
- O 100×50 é grosso demais: perde um atrator.

---

### 2.5 Correção transversal: atratores duplicados e integração numérica (afeta M1 e M2 do SCM e do GCM)

Esta é a parte mais importante da fase, porque mudou o entendimento de um problema que vinha desde o SCM.

**Sintoma.**
- Teste 6 v2: referência HD 1000×1000 com 5 atratores, sendo que o HD tem 4 `[term:10/10]`.
- No relatório parcial, o SCM já tinha mostrado um atrator duplicado no HD, atribuído a "variações numéricas entre arquiteturas de hardware" `[ctx:REL]`.
- No notebook do aluno, a duplicata aparecia já no 200×200; no PC de mesa, nunca no 200×222 `[conv]`.

**Hipóteses, na ordem em que foram testadas:**

| # | Hipótese | Teste | Resultado | Situação |
|---|---|---|---|---|
| H1 | É um defeito só do GCM com m = 1 | `[diag:ref]`: SCM no mesmo grid 1000×1000 | O SCM deu **6** atratores (pares duplicados a 0,0003 e 0,0004 de distância entre centróides); o GCM m = 1 deu 5 (par a 0,0005) | **Refutada:** o SCM também duplica |
| H2 | Atração lenta (multiplicador de Floquet ≈ ±1), um limite intrínseco do SCM | `[diag:tol]`: Jacobiano do mapa de 1 período por diferenças finitas | Os 4 atratores são pontos fixos de período 1, com \|λ\| = 0,7645 nos quatro (= e^(−aT/2), teórico) | **Refutada** |
| H3 | Erro de integração (reltol 1e-3) da ordem do tamanho da célula | `[diag:tol]`: mesmo SCM com Vern9 1e-8 e 1e-10 | 1000×1000: padrão **6** → Vern9 **4** atratores. 200×222: **62,10%** das células mudam de destino (padrão × 1e-8), e `compare_basins` dá 0,8306. 1e-8 × 1e-10: só 0,0541% (200×222) e 0,1892% (1000×1000) mudam, ou seja, o mapa preciso está convergido | **Parcialmente confirmada:** a tolerância piorava muito, mas não era a causa raiz (ver H4) |
| H4 | Discretização do SCM: reiniciar cada passo no centro da célula cria "soluções periódicas agrupadas" em torno de um único atrator (Hsu §8.2) | `[img:aluno-2]`: o aluno rodou o SCM (já com Vern9) em 1000×1111 → 5 atratores (espiral ciano). `[sp:s1]`: posição do ponto fixo real dentro da célula + trajetória contínua | ver abaixo | **Confirmada** |

**Evidência de H4** `[sp:s1]`, HD 1000×1111 com Vern9:

- O SCM encontrou 5 atratores. O #5 é a célula 580757, **imediatamente abaixo** da célula 581757, que contém o ponto fixo real x* = (0,841355; 0,068955).
- x* fica na posição relativa (0,057; 0,036) da sua célula, ou seja, **a 3,6% da altura da célula acima da fronteira** com a célula de baixo.
- O centro da célula de baixo, depois de 11 períodos (10 de transiente + 1), cai numa posição que ainda pertence a ela mesma. Por isso ela vira um ponto fixo espúrio.
- A trajetória **contínua** a partir do #5 termina na célula 581757 (rótulo 4). Os dois são o mesmo atrator.
- ⚠ **Interpretação (ordem de grandeza):** com |λ| = 0,7645 por período, 11 períodos contraem a distância ao ponto fixo por 0,7645¹¹ ≈ 0,052 (derivado). A imagem do centro de uma célula vizinha fica a ~5% de célula de x*. Se x* está a menos que isso de uma fronteira, as duas células podem virar pontos fixos.
- Nos grids 1000×1000 e 200×222, nenhum ponto fixo estava tão perto de uma fronteira, e o SCM com Vern9 deu 4 atratores. Por isso a duplicata **depende do grid**.
- O PoliBoA original não sofre com isso porque segue a **trajetória contínua**, sem reiniciar no centro da célula: toda trajetória converge para x*, que está numa única célula. O aluno confirmou o original com 4 atratores no 1000 `[img:aluno-1]`.

**Por que diferia entre o notebook e o PC.** ⚠ **Interpretação:**
- Com as mesmas versões e o mesmo código, sobra o hardware. O Julia compila para a CPU da máquina (o PC é um i5-11400F, alvo `rocketlake` no LLVM 16.0.6), e instruções diferentes (FMA, largura de vetor) mudam o último bit dos cálculos.
- Com reltol 1e-3, uma diferença de último bit pode inverter decisões de aceitar ou rejeitar passos do integrador adaptativo. Os pontos finais passam a diferir na ordem da tolerância.
- No 200×222, x* do atrator amarelo está a ~0,0012 (9% da célula) de uma fronteira `[sp:s1]`, a mesma ordem desse erro. A diferença de máquina basta para criar ou não a célula espúria.
- **Conclusão:** a explicação de arquitetura do relatório parcial estava **certa, mas incompleta**. O hardware era o gatilho; o mecanismo é o de Hsu §8.2, amplificado pela tolerância frouxa.
- Para documentar, sugiro registrar a CPU do notebook com `versioninfo()`.

**Correção implementada (10/10):**
1. **Integração:** `Vern9(); reltol = 1e-8, abstol = 1e-8` em todos os `init()` do SCM e do GCM. O `poliboa_original.jl` **não foi alterado**, por decisão do aluno: é o grupo controle. **Revertido em 11/10** (ver "Decisão final", no fim desta seção).
2. **Agrupamento de soluções periódicas (Hsu §8.2):**
   - SCM: `merge_clustered_attractors!` ao final de `find_attractors_from_scm`.
   - GCM: `merge_clustered_groups`, dentro de `analyze_gcm`.
   - Funções auxiliares em `desenvolvimento_scm/utils.jl`: `continuous_trajectory_cell` (integra a trajetória contínua por `transient_cycles + attractor_tolerance` períodos, com Vern9 1e-8) e `cluster_representatives`.
   - **Critério:** se a trajetória contínua de um atrator termina na bacia de outro, os dois são unidos. As bacias são reatribuídas e os atratores renumerados.
   - **Por que é seguro:** atratores **distintos** nunca são unidos, porque a trajetória contínua de cada um permanece nele. Isso evita o problema da proposta rejeitada em setembro de "juntar atratores próximos", que apagaria atratores vizinhos do SC `[ctx:HIST]`.

**Validação da correção:**

| Caso | Antes | Depois | Fonte |
|---|---|---|---|
| SCM HD 1000×1111 | 5 atratores (`[img:aluno-2]`, `[sp:s1]`) | 4 (0,2088 / 0,1226 / 0,1222 / 0,5464) | `[sp:s3]` |
| SCM HD 1000×1111, versão serial (`build_scmap`, como no teste da M1) | — | 4, mesmas frações; 107,0 s | `[sp:s5]` |
| SCM HD 1000×1000 | 6 (tol. padrão), `[diag:ref]` | 4 (0,2090 / 0,1233 / 0,5458 / 0,1219) | `[sp:s3]` |
| SCM HD 200×222 | 4 | 4 | `[sp:s3]` |
| SCM SC 200×200 | 2 | 2 | `[sp:s3]` |
| GCM m = 1, HD 1000×1111 | 5 | 4; `compare_basins` com o SCM = 0,9998 | `[sp:s2]`, `[sp:s3]` |
| GCM m = 10, HD 1000×1111 | 4 (as amostras já ligavam as células; grupo de 2 células); 260,2 s | — | `[sp:s2]` |
| Original HD 1000×1111 (repositório, sem alteração) | 4 (0,1209 / 0,5510 / 0,1205 / 0,2076), 37,6 s | — | `[sp:s3]` |
| SCM × Original (padrão), HD 1000×1111 | — | 0,8336 | `[sp:s3]` |
| SCM × Original (padrão), HD 200×222 | — | 0,8291 | `[sp:s3]` |
| Original padrão × **cópia** do original com Vern9 1e-8, HD 1000×1111 | — | **0,8338** | `[sp:s4]` |
| SCM × cópia do original com Vern9, HD 1000×1111 | — | **0,9828** | `[sp:s4]` |
| SCM × cópia do original com Vern9, HD 200×222 | — | **0,9735** | `[sp:s4]` |
| Teste 6, referência HD 1000×1000 | 5 atratores | 4 | `[log:m5]` |

⚠ **Interpretação:** a diferença de ~17% entre SCM e original vinha **inteira** da imprecisão do original. O próprio original muda 16,6% das células quando integrado com precisão, e o SCM corrigido concorda 98,3% com essa versão precisa.

**Custo.** O Vern9 1e-8 deixou o SCM e o GCM cerca de 2,8 a 3,0× mais lentos (derivado):
- GCM HD: 14,88 → 41,90 s;
- GCM SC: 86,42 → 249,93 s;
- SCM HD: 1,43 → 4,27 s;
- SCM SC: 8,48 → 24,55 s.

Fontes: `[term:10/10]`, `[log:m5]`. O agrupamento custa uma integração de ~90 períodos (HD) ou ~180 (SC) por atrator, desprezível.

#### Decisão final (11/10, noite): integrador igual ao do original

**Decisão do aluno.** SCM e GCM devem ficar o mais parecidos possível com o PoliBoA original e recuperar a velocidade. Por isso:

- **O integrador volta a ser o do original.** Todos os `init()` usam exatamente a chamada do `poliboa_original.jl`: `init(ode_problem; dense=false, save_everystep=false, save_start=false, maxiters=1e10)`. O `continuous_trajectory_cell` também passou a usar essa chamada.
- **O agrupamento por trajetória contínua fica.** Ele é apresentado como a **análise de agrupamento das soluções periódicas, parte do método SCM** (Hsu §8.2), e não como correção de erro numérico. Os comentários de `milestone2.jl` foram reescritos nesse sentido.
- **O `reset_dt` do `map_point` (GCM) fica**, porque garante que 1 e 4 threads deem o mesmo mapa, que é a comparação da M4.

**M5 enxuta** (`tests_milestone5_gcm.jl`, `milestone5_gcm.jl`). Agora segue o `desenvolvimento_scm/tests_milestone5.jl`: mesmos testes, mesma ordem (HD e depois SC, Testes 1–5) e mesmos parâmetros (divs = 200, trials = 3, divs_range = [100, 200, 400]).

- **Acréscimos do GCM que ficaram:**
  - Teste 1 com `benchmark_method_gcm` (integrações N·m);
  - Teste 2 com GCM × SCM e concordância por domicílio;
  - Teste 3 com GCM × SCM (custo ≈ m);
  - Teste 4 com Ī e Mapa (MB);
  - Teste 5 com `thread_scaling_gcm`;
  - Teste 6 (refinamento, pedido do orientador) por último.
- **O que saiu:**
  - Teste 0 e `reproducibility_gcm`, porque a M4 já compara 1 × 4 threads;
  - o Teste 6 estratificado (a opção `sampling` continua no `GcmProblem`);
  - a linha m = 200 do SC.
- **O que mudou no Teste 6:**
  - a referência fina é calculada **uma vez por sistema** (`refinement_reference_gcm`) e passada pronta para `refinement_study_gcm`;
  - no HD, cada grid grosso ganha a linha m = K (`full_sampling = true`), que mostra o teto do "ponto a ponto".
- O cabeçalho e o "Resumo Final" foram ajustados: saiu a linha de reprodutibilidade.

⚠ **Correção de interpretação:** a frase "o Vern9 removeu a causa" (H3) deve ser lida assim: o Vern9 **reduzia** o ruído do mapa, mas a causa dos atratores duplicados é a de H4, e quem a resolve é o agrupamento. A revalidação abaixo confirma isso: sem Vern9, o agrupamento sozinho dá o número correto de atratores em todos os casos.

**Revalidação com o integrador do original** (`s6_revalidate.jl`, 4 threads salvo indicação; log `SCRATCH/s6.log`, `[log:s6]`):

| Caso | Atratores | Tempo | Observação |
|---|---|---|---|
| SCM HD 1000×1111 | **4** (0,1240 / 0,5535 / 0,2081 / 0,1144) | 10,7 s | antes do agrupamento eram 5 (`[img:aluno-2]`) |
| SCM HD 1000×1000 | **4** (0,5528 / 0,1239 / 0,1153 / 0,2080) | 9,9 s | antes eram 6 (`[diag:ref]`) |
| SCM HD 200×222 | **4** | 0,4 s | |
| SCM SC 200×200 | **2** (−1: 0,0266) | 2,7 s | |
| GCM m = 1, HD 1000×1111 | **4** | 11,1 s | `compare_basins` com o SCM = 0,8667 (ver §8) |
| GCM m = 10, HD 200×222, 1 thread | **4** | **14,8 s** | referência de 10/10: ≈ 15 s |
| GCM m = 10, SC 200×200, 1 thread | **2** | **83,8 s** | referência de 10/10: ≈ 86 s |
| Original HD 1000×1111, 1 thread | 4 | 31,8 s | SCM × original = **0,8715** |
| Original HD 200×222, 1 thread | 4 | 1,2 s | SCM × original = **0,8678** |

As imagens do SCM 1000×1111 e 1000×1000 e do GCM SC foram conferidas visualmente: 4 bacias no HD e 2 no SC. A região externa do SC fica granulada, como no original (`SCRATCH/s6_*.png`).

**Testes oficiais com a versão final** (11/10, noite; 4 threads; `[log:m1b]`, `[log:m4b]`, `[log:m5b]` = `SCRATCH/m1b.log`, `m4b.log`, `m5b.log`):

- **M1** (132 s):
  - HD com **4** atratores: 4947 / 5425 / 9332 / 24696;
  - SC com **2**: 15495 / 24212; 293 divergentes (derivado).
  - Esses números são idênticos aos da v2 de 10/10: no 200×222 não havia duplicata, e o agrupamento não mudou nada.
- **M4** (166 s): imagens da M1 e da M4 **idênticas byte a byte** (HD e SC).
- **M5 enxuta** (2456 s ≈ 41 min): resultados na tabela abaixo.

**Comparação com a rodada Vern9 de 11/10 de manhã** (`[log:m5]` × `[log:m5b]`), só para os testes que existem nas duas versões:

| Teste | Métrica | Vern9 1e-8 | Integrador do original |
|---|---|---|---|
| 1 | Tempo médio do GCM, HD / SC | 41,90 s / 249,93 s | **14,94 s / 84,15 s** |
| 2 HD | GCM × Original; GCM × SCM | 83,43%; 86,76% | 84,73%; 85,08% |
| 2 HD | GCM × SCM, domicílio único; fronteira | 99,97%; 35,39% | 99,75%; 36,17% |
| 2 SC | GCM × Original; GCM × SCM | 22,85%; 95,71% | 26,91%; **74,92%** |
| 2 SC | GCM × SCM, domicílio único; fronteira | 99,47%; 27,23% | 96,62%; **60,76%** |
| 3 | Tempo GCM/SCM, HD / SC | 10,06× / 10,21× | 10,25× / 9,86× |
| 3 | Tempo GCM/Original, HD / SC | 31,82× / 29,85× | **11,93× / 10,26×** |
| 4 | MB/kCell HD (100 / 200 / 400) | 11,813 / 11,728 / 11,663 | 11,816 / 11,742 / 11,682 |
| 4 | Ī SC (100 / 200 / 400) | 2,618 / 2,976 / 2,998 | 3,547 / 4,603 / 5,226 |
| 5 | Speedup com 4 threads, HD / SC | 3,39× / 3,26× | 3,19× / 3,13× |
| 6 HD | Referência 1000×1000 | 4 atratores (0,5458 / 0,2090 / 0,1233 / 0,1219); 27,5 s | **4** (0,5528 / 0,2073 / 0,1237 / 0,1162); 11,3 s |
| 6 HD | 200×100, m = 5: erro fração; erro/célula; ponto a ponto | 0,38 p.p.; 0,0948; 81,01% | 0,24 p.p.; 0,0960; 80,56% |
| 6 HD | 200×100, m = 10 | 0,35 p.p.; 0,0631; 82,22% | 0,14 p.p.; 0,0643; 81,51% |
| 6 HD | 200×100, m = 20 | 0,32 p.p.; 0,0402; 82,77% | 0,09 p.p.; 0,0406; 82,21% |
| 6 HD | 200×100, m = K = 50 | 0,21 p.p.; 0,0094; 83,22% | 0,04 p.p.; 0,0103; 82,60% |
| 6 HD | 100×50, m = 5 / 10 / 20 / 200 | 4 de 4 / 3 / 3 / 3 de 4 | 3 de 4 em todos (erro 11,62 p.p. = bacia de 0,1162 perdida) |
| 6 SC | Referência 1000×1000 | 2 (−1: 0,0260); 188,7 s | 2 (−1: 0,0257); 66,0 s |
| 6 SC | 200×100, m = 5 | 0,04 p.p.; 0,0670; 87,81% | 0,26 p.p.; 0,1194; 73,89% |
| 6 SC | 100×50, m = 5 | 0,11 p.p.; 0,1019; 80,68% | 0,22 p.p.; 0,1456; 68,20% |

Como ler a comparação:
- **Correção.** A referência HD do Teste 6 dava 5 atratores com o integrador do original em 10/10 (`[term:10/10]`). Agora dá **4**: a bacia de 0,1162 é a soma das duas bacias duplicadas da rodada de 10/10 (0,0244 + 0,0918). Isso mostra que o agrupamento corrige o problema sem precisar do Vern9.
- **Velocidade.** Os tempos voltaram aos de 10/10 (Teste 1: HD 14,94 s; SC 84,15 s).
- ⚠ **Interpretação (SC).** Com o integrador do original, o SC fica bem mais "ruidoso":
  - a fronteira sobe para 60,76% das células;
  - GCM × SCM cai para 74,92%;
  - Ī cresce com N;
  - as métricas do Teste 6 pioram (ponto a ponto de 73,89% contra 87,81%).
  
  É o mesmo efeito documentado no Anexo A: com 100 períodos de transiente e reltol 1e-3, amostras vizinhas caem em destinos diferentes por erro numérico. Não é erro do método; é a precisão do integrador escolhido, a mesma do original.

---

# Parte 2 — Diagnósticos, desafios e alterações no SCM

## 3. Diagnósticos e testes

| ID | Objetivo | Configuração | Resultado | Conclusão | Fonte |
|---|---|---|---|---|---|
| P1 | Experimento do orientador (protótipo) | Python, RK4 de 200 passos por período; HD com fino 1000×1000; grossos 100×50, 200×100 e 500×500 | Referência com 4 atratores; 200×100 m = 5: TV 0,0954, ponto a ponto 81,18%; 100×50 perde 1 atrator | Base dos números esperados | `[ctx:PROTO]` |
| P2 | Estratificada × aleatória | `strat_test.py` | ≈ 13% menos erro por célula | Motivou `sampling = :stratified` | `[ctx:CLAUDE]` |
| P3 | Rodadas adaptativas | `adapt_test.py`: 5 amostras em todas as células + 50 só nas de domicílio múltiplo | TV 0,021 com 434k integrações, contra 0,040 com 400k (m = 20 uniforme) | Ideia promissora, **não implementada** em Julia | `[ctx:CLAUDE]` |
| V1-a | M1 v1 em Julia (setembro) | HD e SC 200, S = 3 | HD 4, SC 2 | OK | `[ctx:CLAUDE]` |
| V1-b | M4 v1 | 1 × 4 threads | Diferentes (mais células no verde) | Levou ao `reset_dt` | `[ctx:CLAUDE]` |
| D1 | Repositório × v2 do pacote | `diff` ignorando fim de linha | Conteúdo idêntico; só CRLF e quebra final | v2 já estava aplicada | `[conv]` |
| D2 | M1 v2 | `tests_milestone1_gcm.jl` | HD 4; SC 2 (293 divergentes) | OK | `[term:10/10]` |
| D3 | M4 v2 | `tests_milestone4_gcm.jl` + `cmp` das imagens | M1 = M4 byte a byte | Reprodutível | `[conv]` |
| D4 | M5 v2 completa | `tests_milestone5_gcm.jl` | Teste 0 OK; referência do Teste 6 HD com **5** atratores; SC × Original 26,9% | Abriu duas investigações | `[term:10/10]` |
| D5 | Revisão de código (M2, M4, M5, testes) | Leitura contra Hsu e contra o SCM | Lógica correta; 5 pontos (1 cosmético corrigido) | — | `[conv]` |
| D6 | Origem do 5º atrator; SC × Original | `diagnostic_ref_gcm.jl` | GCM m = 1: 5 atratores; SCM: 6; `compare_basins` = 0,8114. SC: SCM × Original = 0,2726; divergentes: SCM 2,50%, Original 62,92%, original com 1 atrator | Não é defeito só do GCM; original falha no SC | `[diag:ref]` |
| D7 | Tolerância + multiplicadores de Floquet | `diagnostic_tol_gcm.jl` | \|λ\| = 0,7645 nos 4; 1000×1000: 6 → 4 com Vern9; mapa convergido entre 1e-8 e 1e-10 | Tolerância piora muito; atração lenta descartada | `[diag:tol]` |
| D8 | Original × SCM em 1000 (aluno) | `poliboa.jl` com `test_helmduff(1000)`; `tests_milestone1.jl` com `test_scm_helmduff(1000)` (já com Vern9) | Original: 4; SCM: 5 (espiral ciano) | Tolerância não era a causa raiz | `[img:aluno-1]`, `[img:aluno-2]` |
| D9 | Posição de x* na célula; trajetória contínua | `s1_scm_1000.jl` (1000×1111, 1000×1000, 200×222) | #5 = célula vizinha de x*; x* a 3,6% da fronteira; trajetória contínua do #5 termina no #4 | **Causa raiz confirmada** (Hsu §8.2) | `[sp:s1]` |
| D10 | GCM antes da correção, 1000×1111 | `s2_gcm_before.jl`, m = 1 e m = 10 | m = 1: 5 atratores; m = 10: 4 (260,2 s) | GCM com m > 1 se autocorrige | `[sp:s2]` |
| D11 | Validação da correção + original | `s3_after_fix.jl` | SCM 4/4/4 (HD), 2 (SC); GCM m = 1: 4; SCM × Original: 0,83 | Correção funciona | `[sp:s3]` |
| D12 | Quanto da diferença SCM × original vem do original | `s4_original_tol.jl` (**cópia** do original com Vern9) | Original padrão × Vern9: 0,8338; SCM × original Vern9: 0,9828 (1000×1111) e 0,9735 (200×222); SC: original Vern9 ainda com 1 atrator | A diferença vem do original | `[sp:s4]` |
| D13 | Região do SC; equivalentes M1, Teste 0 e Teste 6; SCM serial | `s5_validation.jl` | SC [−π, π]: original sem alteração com 2 atratores e 0% divergentes; GCM M1: HD 4, SC 2; reprodutível; Teste 6 HD OK; SCM serial 1000×1111: 4 | Região [−3,14; 3,14] causa a falha do original no SC | `[log:s5]` (`SCRATCH/s5.log`) |
| D14 | M1 oficial, v2.1 | `tests_milestone1_gcm.jl` (aluno, 11/10) | HD 4; SC 2 (466 divergentes) | OK | `[term:11/10]` |
| D15 | M4 oficial, v2.1 | `tests_milestone4_gcm.jl` (Claude, 11/10) | 7 min 14 s; imagens M1 = M4 | OK | `[log:m4]` |
| D16 | M5 oficial, v2.1 | `tests_milestone5_gcm.jl` (Claude, 11/10) | 128 min 54 s; todos os resultados da §2.4 | OK | `[log:m5]` |

---

## 4. Desafios encontrados e explorações feitas

### 4.1 O que é um "GCM correto" (v0 → v1)

- A primeira M2 (v0) usava só a imagem de maior probabilidade, o que equivale ao "SCM dominante" de Hsu §11.3.
- Ela foi descartada porque não usa as probabilidades, que são a razão de existir do GCM `[ctx:HIST]`.
- A v1 implementou a teoria completa: grupos persistentes, período, domicílios e absorção.
- **O que não funcionou:** a v0 da M2.

### 4.2 Amostragem determinística × sorteio (v1 → v2)

- A v1 usava um retículo fixo S×S (S = 3).
- O orientador propôs dividir a célula em K subcélulas e sortear m sem reposição, por rodadas `[ctx:CLAUDE]`. Isso permite:
  - (a) variar o custo continuamente via m;
  - (b) alinhar as amostras ao grid fino, que é a base do Teste 6.
- **O que foi substituído:** o retículo S×S. Ele não falhou; foi substituído por uma estratégia mais flexível.

### 4.3 Dependência do número de threads (v1)

- Hipótese: o passo adaptativo do integrador era herdado entre integrações.
- Correção: `reset_dt = true` + semente por célula.
- **Confirmada** pelo Teste 0 (0 diferenças) e pela igualdade byte a byte das imagens.

### 4.4 Atratores duplicados (§2.5)

**O que não funcionou ou foi abandonado:**
- "Sensibilidade da máquina" como explicação completa (relatório parcial): incompleta.
- "Juntar atratores próximos" (setembro): rejeitada, porque apagaria atratores próximos do SC `[ctx:HIST]`.
- Trocar só a tolerância (D7): reduziu o problema (6 → 4 no 1000×1000), mas não eliminou (5 no 1000×1111, D8).

**O que funcionou:** o agrupamento por trajetória contínua (Hsu §8.2).

### 4.5 Confusão de métricas

- "Ponto a ponto" do Teste 6 (~81%) × concordância no mesmo grid (~98%): ver o quadro de atenção na §2.4.
- "Bytes alocados" (MB/kCell ≈ 11,7) × memória do mapa (MB): ver o Teste 4.
- O 86,76% de GCM × SCM no Teste 2 v2.1 coincide numericamente com o 86,8% de SCM × original do relatório parcial, mas **são comparações diferentes**.

### 4.6 Spring Column × original

- O original dava 1 atrator e ~63% de células divergentes no SC. Isso acontecia com reltol 1e-3 (62,92%, D6) e também com Vern9, na cópia (64,87%, D12).
- O relatório parcial atribuía isso a o original marcar a trajetória inteira como divergente. O mecanismo existe, mas **o gatilho** é a região: x₁ ∈ [−3,14; 3,14], e não [−π; π].
  - Sobra uma faixa de 0,0016 em cada ponta, fora da região principal **e** da estendida (que é a mesma em x₁).
  - Um ponto que cai nela é tratado como divergente.
  - ⚠ **Interpretação do mecanismo:** como o atrator em x₁ ≈ 3,12 fica perto dessa ponta, as trajetórias do original que se aproximam dele passam pela faixa e são inteiramente marcadas como divergentes.
- **Evidência (D13):** com [−π; π], o original **sem nenhuma alteração** encontra 2 atratores e 0% de divergentes; o SCM também.
- **Decisão do aluno (11/10):** manter [−3,14; 3,14], igual ao original, e explicar a faixa em rodapé no relatório.
- **Consequência:** SCM e GCM mostram ~2,6% e ~1,2% de células "divergentes" no SC, que são artefato da faixa. É a mesma origem das manchas pretas em `gcm_basin_sc_testm{1,4}.png`.

### 4.7 Regras de trabalho que surgiram de erros

Do `[ctx:HIST]`:
- testes só trocam nomes (o assistente tinha inventado testes);
- não mudar algoritmo sem avisar;
- comparações entre métodos só na M5.

Desta fase:
- o original não pode ser alterado (grupo controle);
- toda conclusão precisa de evidência antes de ser afirmada. A hipótese da tolerância (H3) chegou a ser apresentada como causa antes de ser testada contra o grid 1000×1111, e o aluno corrigiu.

---

## 5. Alterações no SCM

| Arquivo | O que mudou | Por quê | Quando |
|---|---|---|---|
| `milestone5.jl` | `compare_basins`: atribuição gulosa **bijetiva**, com −1 ↔ −1 fixo e `get(mapping, x, nothing)` | A versão antiga casava −1 com o atrator amarelo no SC | `[git:15f367b]`, 15/08 |
| `utils.jl`, `milestone1.jl`, `milestone4.jl` | `compute_single_cell_mapping` movido para `utils.jl` (desduplicação) | Havia cópias divergentes | `[git:15f367b]` |
| `milestone1.jl` (`build_scmap`) e `milestone4.jl` (`process_cell_range`, `build_scmap_parallel`) | Vern9 1e-8 em 10/10; **revertido em 11/10** para a chamada exata do original. Hoje estão **idênticos ao commit `fdc51bc`** (`git diff` vazio) | Decisão do aluno: mesmo integrador do original e velocidade | 10/10 → 11/10 |
| `utils.jl` | **Novas:** `continuous_trajectory_cell(u0, bp)`, que agora usa o mesmo `init()` do original + `step!` (sem Vern9), e `cluster_representatives(target)` | Análise de agrupamento das soluções periódicas (Hsu §8.2) | 10/10 (ajustado em 11/10), não commitado |
| `milestone2.jl` | **Nova:** `merge_clustered_attractors!(grid, attractors, bp)`, chamada no fim de `find_attractors_from_scm`. Comentários reescritos em 11/10 para apresentar a etapa como parte do método SCM (Hsu §8.2) | Idem | 10/10 (ajustado em 11/10), não commitado |
| `poliboa_original.jl` | **Nenhuma alteração** | Grupo controle (decisão do aluno) | — |

**Impacto no comportamento (versão final):**

- **A única diferença em relação ao SCM do relatório parcial é o agrupamento.** O integrador é o mesmo. Nos grids em que não há soluções periódicas duplicadas (HD 200×222, SC 200×200), o resultado é igual ao de antes: a M1 de 11/10 à noite reproduziu os números de 10/10.
- **Em grids finos, o agrupamento muda o número de atratores para o correto:** HD 1000×1111 de 5 para 4; HD 1000×1000 de 6 para 4 (`[log:s6]`).
- **SCM × original:** 86,78% no 200×222 e 87,15% no 1000×1111 (`[log:s6]`), na faixa dos 86,8% do relatório parcial.
- **Tempo:** sem mudança relevante. O agrupamento custa uma integração contínua por atrator.
- **A numeração dos atratores** pode mudar quando há agrupamento: os atratores restantes são renumerados 1, 2, 3…
- ⚠ A M5 **do SCM** (`desenvolvimento_scm/tests_milestone5.jl`) **não foi rodada** nesta fase. Como o integrador voltou a ser o mesmo, os números do relatório parcial devem se manter no 200×222. Vale rodar para confirmar.

**Retrocompatibilidade:**

- Nenhuma assinatura pública mudou: `build_scmap`, `build_scmap_parallel`, `find_attractors_from_scm` e `compute_single_cell_mapping` têm os mesmos argumentos e retornos. O `BasinResult` é o mesmo.
- Os testes do SCM rodam sem alteração. O aluno rodou `tests_milestone1.jl`, e o `s5` rodou o equivalente.
- O `poliboa_original.jl` usa funções de `utils.jl` (`store_cell_center!`, `set_integrator!`, `get_cell_number` etc.), que **não foram alteradas**; só foram acrescentadas funções novas. Por isso o comportamento do original é o mesmo, e as execuções do original no repositório reproduziram os 62,92% de divergentes no SC (D6 e D13).
- Não há requisito novo de pacote: o Vern9 saiu do código.

**O que mais poderia ficar diferente do original** (sugestões; nada disso foi alterado):

1. **`reset_dt` também no SCM.** O `compute_single_cell_mapping` do SCM reaproveita o passo adaptativo da célula anterior. Com o integrador do original, isso pesa:
   - no HD 1000×1111, o GCM com m = 1 (mesmo mapa, mas com `reset_dt`) concorda só **86,67%** com o SCM (`[log:s6]`);
   - é provável que a M4 do SCM (1 × 4 threads) não dê resultados idênticos.
   
   Colocar `reset_dt` no SCM seria a mesma correção que o GCM já tem. Não mexe no original, que segue trajetórias contínuas.
2. **Região do SC [−π; π].** Eliminaria os ~2,6% de células "divergentes" artificiais (§4.6). A decisão foi manter [−3,14; 3,14], igual ao original, com rodapé.
3. **Tolerância do integrador.** Mantida igual à do original por decisão. O custo em precisão está registrado no Anexo A.

---

# Parte 3 — Inventário, reprodutibilidade, estado e anexos

## 6. Inventário de arquivos

Situação em relação ao último commit (`fdc51bc`), segundo `git status` em 11/10/2026.

### 6.1 GCM (`desenvolvimento_gcm/`)

| Arquivo | Situação no git | Finalidade | Origem do conteúdo atual |
|---|---|---|---|
| `structures_gcm.jl` | alterado | `GcmProblem`, `GeneralizedCellMap` (CSR), `PreImageArray`, `GcmAnalysis` | v2 |
| `utils_gcm.jl` | alterado | sorteio, amostras, `map_point` (com `reset_dt`), integrador, acessores, relatórios | v2; o Vern9 entrou em 10/10 e saiu em 11/10, e `make_gcm_integrator` usa a chamada do original |
| `milestone1_gcm.jl` | alterado | `build_gcmap`, `compute_cell_transition_probabilities` | v2 |
| `milestone2_gcm.jl` | alterado | análise de Markov de Hsu + `merge_clustered_groups` (Hsu §8.2) | v2 + agrupamento |
| `milestone4_gcm.jl` | novo (não versionado) | paralelização | v1/v2 |
| `milestone5_gcm.jl` | novo | validação. **Enxugado em 11/10:** saiu `reproducibility_gcm`; entrou `refinement_reference_gcm`; `refinement_study_gcm` recebe a referência pronta e tem a opção `full_sampling` (linha m = K) | v2 + M5 enxuta |
| `tests_milestone1_gcm.jl` | alterado | teste M1 (só nomes trocados em relação ao SCM) | v2 |
| `tests_milestone4_gcm.jl` | novo | teste M4 (idem) | v2 |
| `tests_milestone5_gcm.jl` | novo | **Enxugado em 11/10:** Testes 1–5 como no SCM (HD, depois SC) + Teste 6 (só aleatório; HD com linha m = K) | v2 + M5 enxuta |
| `gcm_basin_{hd,sc}_testm1.png` | alterados | figuras da M1 | regeneradas em 11/10 à noite (versão final) |
| `gcm_basin_{hd,sc}_testm4.png` | novos | figuras da M4 | regeneradas em 11/10 à noite; idênticas às da M1 |
| ~~`diagnostic_ref_gcm.jl`~~, ~~`diagnostic_tol_gcm.jl`~~ | **apagados em 11/10** | diagnósticos D6 e D7 | nenhum outro arquivo dependia deles; a investigação está resolvida; o estudo de precisão fica só no Anexo A |
| ~~`diag_ref_gcm_hd_1000.png`~~, ~~`diag_ref_scm_hd_1000.png`~~, ~~`diag_tol_scm_hd_1000.png`~~ | **apagados em 11/10** | figuras de D6 e D7 | idem |

### 6.2 SCM (`desenvolvimento_scm/`)

| Arquivo | Situação no git | Finalidade da mudança |
|---|---|---|
| `milestone1.jl` | **igual ao commit** | Vern9 entrou em 10/10 e foi revertido em 11/10 |
| `milestone2.jl` | alterado | `merge_clustered_attractors!` (análise de agrupamento, Hsu §8.2) |
| `milestone4.jl` | **igual ao commit** | Vern9 entrou em 10/10 e foi revertido em 11/10 |
| `utils.jl` | alterado | `continuous_trajectory_cell` (integrador do original), `cluster_representatives` |
| `scm_basin_{hd,sc}_testm1.png` | novos | gerados pelo aluno em 10/10 às 18:07, com `test_scm_helmduff(1000)`. A figura do HD é a de 5 atratores **antes** da correção (`[img:aluno-2]`). ⚠ Rodar o teste de novo para atualizar |
| `poliboa_original.jl`, `structures.jl`, `milestone5.jl`, `tests_milestone*.jl`, `diagnostic.jl` | sem alteração nesta fase | — |

### 6.3 Outros arquivos

| Arquivo | Situação | Observação |
|---|---|---|
| `basin_hd.csv` (raiz) | novo, 2,2 MB, 10/10 às 18:15 | ⚠ Interpretação: provavelmente gerado pelo `ZZ - Poli BoA og/poliboa.jl` do aluno, na execução 1000×1111 (o tamanho bate com ~1,1 milhão de rótulos). Não usado aqui |
| `ZZ - Poli BoA og/poliboa.jl`, `desenvolvimento_scm/tests_milestone1.jl` | sem alteração em relação ao commit | O aluno trocou 200 → 1000 em 10/10 e depois desfez |

### 6.4 Fora do repositório (`SCRATCH`)

- **Scripts:** `common.jl`, `s1_scm_1000.jl`, `s2_gcm_before.jl`, `s3_after_fix.jl`, `s4_original_tol.jl`, `s5_validation.jl`, `s6_revalidate.jl` (revalidação com o integrador do original), `s7_smoke_m5.jl` (teste rápido do Teste 6 novo), `poliboa_original_vern9_COPIA.jl` (cópia **modificada** do original, só para D12).
- **Logs:**
  - rodada Vern9: `s5.log`, `m4.log`, `m5.log`;
  - versão final: `s6.log`, `m1b.log`, `m4b.log`, `m5b.log`.
- **Resultados serializados:** `scm_*.jls`, `gcm_before_*.jls`, `orig_hd_1000x1111.jls`.
- **Figuras:** ver §9.

---

## 7. Reprodutibilidade

**Ambiente (PC de mesa, 11/10/2026):**

| Item | Valor |
|---|---|
| Julia | 1.11.6 (commit 9615af0f26) |
| SO | Windows 11 Home 10.0.26200 (x86_64-w64-mingw32) |
| CPU | 12 × 11th Gen Intel Core i5-11400F @ 2,60 GHz (LLVM 16.0.6, alvo `rocketlake`) |
| DifferentialEquations | v7.16.1 |
| Plots | v1.40.19 |
| ProgressMeter | v1.11.0 |
| Accessors | v0.1.42 |
| Outros pacotes (`CLAUDE.md`) | PlotThemes, LaTeXStrings, Colors, Statistics, Random, SparseArrays |

Fonte: `versioninfo()` e `Pkg.status` executados em 11/10.

⚠ O notebook do aluno não teve `versioninfo()` registrado.

**Comandos** (na pasta `desenvolvimento_gcm/`, PowerShell; o `&&` não existe no PowerShell 5.1):

```bash
julia -t 4 tests_milestone1_gcm.jl
```

```bash
julia -t 4 tests_milestone4_gcm.jl
```

```bash
julia -t 4 tests_milestone5_gcm.jl
```

Tempos medidos com 4 threads:

| Versão | M1 | M4 | M5 | Fonte |
|---|---|---|---|---|
| Final (integrador do original, M5 enxuta) | 132 s | 166 s | 2456 s ≈ 41 min | `[log:m1b]`, `[log:m4b]`, `[log:m5b]` |
| Vern9 (11/10 de manhã, M5 antiga) | não registrado | 7 min 14 s | 128 min 54 s | `[log:m4]`, `[log:m5]` |

**Configurações dos sistemas** (iguais em todos os testes):

- **HD:**
  - parâmetros [0,1, −1,2, −0,3, 2,0, e = 0,077, ω = 1,17];
  - período 2π/1,17; 10 períodos de transiente; `maximum_cycles` = 1000; `maximum_extended_cycles` = 20; `attractor_tolerance` = 80;
  - região [−1,2; 1,5] × [−1,5; 1,5], estendida [−2,0; 2,5] × [−2,0; 2,0];
  - grid padrão 200 × ⌊(3,0/2,7)·200⌋ = 200×222.
- **SC:**
  - parâmetros [0,01, 0,05, 0,8, 0,01, 0,05, ω = 0];
  - período 7,853981634; **100** períodos de transiente;
  - região [−3,14; 3,14] (cíclica) × [−1,3; 1,3], estendida [−3,14; 3,14] × [−2,0; 2,0];
  - grid 200×200.
- **GCM:** subgrid 10×10, m = 10, `sampling = :random`, `seed = 1`. O sorteio usa `Xoshiro(seed·10⁹ + cell_id)`.
- **Integrador (SCM, GCM e original):** o mesmo nos três. É a chamada do `poliboa_original.jl`, `init(ode_problem; dense=false, save_everystep=false, save_start=false, maxiters=1e10)`, que usa o padrão do DifferentialEquations (reltol 1e-3, abstol 1e-6). O GCM reinicia o integrador a cada amostra (`reset_dt = true`).
- **M5:** divs = 200, trials = 3, divs_range = [100, 200, 400], como no SCM. Teste 6 com fino 1000×1000, grossos [[100, 50], [200, 100]], m ∈ {5, 10, 20} (+ m = K no HD), sorteio aleatório.

**Versões do código:**

| Referência | O que contém |
|---|---|
| `main` @ `fdc51bc` | GCM v0 + SCM com `compare_basins` corrigido |
| `origin/gcm-milestones-2-4-5` @ `5e48742` | GCM v1 (não inspecionado localmente) |
| Árvore de trabalho local (não commitada) | GCM v2.1 + alterações do SCM desta fase |

⚠ **Recomendação:** commitar o estado atual, para que os números deste documento tenham um commit de referência.

---

## 8. Estado atual, limitações e próximos passos

**Estado (versão final, 11/10 à noite).** O GCM está completo nas quatro milestones, com o **mesmo integrador do PoliBoA original**:
- HD com 4 atratores e SC com 2;
- M1 = M4 byte a byte;
- memória O(N);
- Teste 6 sem atrator duplicado na referência.

O SCM difere do relatório parcial só pela **análise de agrupamento (Hsu §8.2)**, que agora faz parte do método e dá o número correto de atratores também em grids finos. O original foi preservado como grupo controle.

**Limitações conhecidas:**

1. **Custo:** o GCM custa ~10× o SCM com m = 10 (tempo 10,25× no HD e 9,86× no SC, `[log:m5b]`).
2. **Precisão do integrador (o mesmo do original):**
   - mapas ruidosos no SC: fronteira de 60,76%, GCM × SCM de 74,92%, Ī crescendo com N (3,5 → 5,2);
   - sensibilidade ao estado do integrador e à máquina.
   
   Ver o Anexo A para o que muda com integração precisa.
3. **Teste 6, grid 100×50:** grosso demais para o HD; perde o atrator de bacia 0,1162 em todos os m.
4. **SC:** a faixa morta da região [−3,14; 3,14] gera ~2,7% (SCM) e ~0,7% (GCM) de células "divergentes" artificiais (`[log:s6]`). A decisão é manter e explicar em rodapé.
5. **Critério de agrupamento:** validado só com atratores de período 1 (HD e SC), agora com os dois integradores. ⚠ Não foi testado com atratores caóticos ou quase-periódicos.
6. **Empate exato no argmax** (raro) manda a célula para o sumidouro.
7. ⚠ **`compute_single_cell_mapping` do SCM continua sem `reset_dt`.** Com o integrador do original, o reinício do integrador muda 13,3% dos rótulos no HD 1000×1111 (GCM com m = 1 × SCM = 0,8667, `[log:s6]`). A reprodutibilidade 1 × N threads do **SCM** não foi medida. Sugestão na §5.
8. As rodadas adaptativas (P3) não foram implementadas.

**Próximos passos sugeridos:**

1. Commitar o estado atual (git é com o aluno).
2. Rodar a M5 do **SCM** (`desenvolvimento_scm/tests_milestone5.jl`) para confirmar que os números do relatório parcial se mantêm (o integrador é o mesmo; só entrou o agrupamento).
3. Decidir sobre o `reset_dt` no SCM (§5, sugestão 1).
4. Rodar `tests_milestone1.jl` do SCM com 200 no notebook e registrar `versioninfo()` das duas máquinas.
5. Levar ao orientador a tabela custo × precisão do Teste 6 (§2.5, "Decisão final"), com a ressalva sobre o "ponto a ponto" e, se útil, o Anexo A.
6. Avaliar a implementação das rodadas adaptativas (P3).
7. Relatório final:
   - incluir a análise de agrupamento (Hsu §8.2) como parte do método;
   - incluir o rodapé sobre a faixa da região do SC;
   - substituir a explicação de "arquitetura" pela explicação completa.

---

## 9. Figuras e tabelas

### 9.1 Figuras existentes

| Figura | Caminho | Conteúdo |
|---|---|---|
| Bacias GCM HD 200×222 (M1/M4, v2.1) | `desenvolvimento_gcm/gcm_basin_hd_testm1.png` (= `…testm4.png`) | 4 atratores |
| Bacias GCM SC 200×200 (M1/M4, v2.1) | `desenvolvimento_gcm/gcm_basin_sc_testm1.png` (= `…testm4.png`) | 2 atratores; manchas pretas = faixa morta da região |
| Antes da correção, HD 1000×1000 | `desenvolvimento_gcm/diag_ref_gcm_hd_1000.png`, `diag_ref_scm_hd_1000.png` | GCM m = 1 com 5 atratores; SCM com 6 |
| SCM HD 1000×1000 com Vern9 1e-10 | `desenvolvimento_gcm/diag_tol_scm_hd_1000.png` | 4 atratores |
| Original HD 1000 (aluno) | `[img:aluno-1]` | 4 atratores |
| SCM HD 1000×1111 antes da correção (aluno) | `[img:aluno-2]`; `desenvolvimento_scm/scm_basin_hd_testm1.png` | 5 atratores (espiral ciano) |
| GCM 1000×1111 antes da correção | `SCRATCH/gcm_before_hd_1000_m1.png`, `…_m10.png` | m = 1: 5; m = 10: 4 |
| SCM corrigido | `SCRATCH/scm_fixed_hd_1000x1111.png`, `…_1000x1000.png`, `…_200x222.png`, `scm_fixed_sc_200x200.png`, `scm_serial_hd_1000x1111.png` | 4 / 4 / 4 / 2 / 4 |
| GCM m = 1 corrigido | `SCRATCH/gcm_m1_fixed_hd_1000x1111.png` | 4 |
| Original HD 1000×1111 | `SCRATCH/orig_hd_1000x1111.png` | 4 (mais "granulado" por causa da reltol 1e-3) |
| SC com região [−π; π] | `SCRATCH/scm_fixed_sc_200_pi.png`, `orig_sc_200_pi.png` | 2 atratores, sem divergentes; o original sai granulado |

### 9.2 Tabelas prontas neste documento

- Testes 1–6 da M5 (§2.4).
- Hipóteses dos atratores duplicados (§2.5).
- Validação da correção (§2.5).
- Diagnósticos (§3).

### 9.3 Sugestões de figuras para o relatório

1. **Fluxograma do GCM:** amostragem → mapa CSR → pré-imagens → SCM dominante → grupos persistentes → agrupamento (§8.2) → domicílios → absorção → bacias. Na mesma linha do fluxograma das 4 etapas que o orientador pediu para o relatório parcial `[ctx:HIST]`.
2. **Antes × depois da correção:** SCM 1000×1111 com 5 atratores × corrigido com 4 × original, lado a lado.
3. **Esquema do mecanismo de Hsu §8.2:** duas células vizinhas, x* a 3,6% da fronteira e as imagens dos centros após 11 períodos.
4. **Curva custo × precisão do Teste 6:** erro/célula × integrações, aleatória × estratificada, para 200×100 e 100×50, HD e SC.
5. **Mapa de domicílios (HD):** células de domicílio único × múltiplo, mostrando que a fronteira coincide com as regiões fractais.
6. **Speedup e memória:** barras do Teste 5 e retas MB/kCell × N do Teste 4.
7. **SC, região [−3,14; 3,14] × [−π; π]:** para o rodapé.

---

## 10. Glossário

| Sigla / termo | Significado |
|---|---|
| **SCM** | Simple Cell Mapping: cada célula tem uma única imagem |
| **GCM** | Generalized Cell Mapping: cada célula tem várias imagens, com probabilidades (cadeia de Markov) |
| **CCM** | Cell-to-Cell Mapping, família de métodos de Hsu (1987) |
| **HD** | Oscilador de Helmholtz-Duffing (sistema de teste, 4 atratores) |
| **SC** | Spring Column, coluna com mola (sistema de teste, 2 atratores) |
| **M1, M2, M4, M5** | Milestones: construção do mapa, detecção, paralelização, validação (a M3 foi descartada) |
| **CSR** | Compressed Sparse Row: listas concatenadas + vetor de deslocamentos (`offsets`) |
| **Célula sumidouro** | Célula fictícia (nó N+1) que recebe as trajetórias que saem da região; aparece como −1 nas bacias |
| **SCM dominante** | SCM obtido do GCM tomando, em cada célula, a imagem de maior probabilidade (Hsu §11.3) |
| **Grupo persistente (PG)** | Conjunto fechado de células que se comunicam; é o atrator no GCM |
| **Célula transiente** | Célula que não pertence a nenhum grupo persistente |
| **Período de um grupo** | Número de subgrupos visitados ciclicamente (Hsu §11.6) |
| **Domicílio, Dm(z)** | Grupo persistente em que a célula z pode terminar; domicílio múltiplo = fronteira entre bacias |
| **Probabilidade de absorção** | Probabilidade de a célula terminar em cada grupo; calculada por (I − Q)X = B |
| **Ī** | Número médio de imagens por célula |
| **K, m** | K = número de subcélulas por célula (prod(subgrid)); m = subcélulas sorteadas por célula |
| **TV** | Distância de variação total; o "erro/célula" do Teste 6 |
| **p.p.** | Pontos percentuais |
| **Multiplicador de Floquet (λ)** | Autovalor do Jacobiano do mapa de um período; \|λ\| < 1 indica atração, e quanto mais perto de 1, mais lenta |
| **Vern9** | Método Runge-Kutta de Verner de 9ª ordem (DifferentialEquations.jl) |
| **reltol / abstol** | Tolerâncias relativa e absoluta do controle de passo adaptativo |
| **reset_dt** | Opção de `reinit!` que descarta o passo adaptativo da integração anterior |
| **JIT** | Compilação just-in-time: o Julia gera código para a CPU da máquina |
| **FMA / AVX** | Instruções de CPU (multiplica-e-soma fundido; vetores largos) que alteram o arredondamento no último bit |
| **Grupo controle** | O `poliboa_original.jl`, mantido sem alteração para comparação |
| **IC / FAPESP** | Iniciação Científica / Fundação de Amparo à Pesquisa do Estado de São Paulo |
| **NoDE** | Grupo Non-linear Dynamics for Engineering (Poli-USP) |

---

## Anexo A — Estudo de precisão do integrador (registro; nenhum script ficou no repositório)

Este anexo guarda o que se aprendeu com a integração precisa (Vern9, reltol = abstol = 1e-8 e 1e-10). O Vern9 **não** é usado na versão final, que segue o integrador do original por decisão do aluno.

**Como o estudo foi feito.** Scripts `diagnostic_tol_gcm.jl` (apagado do repositório em 11/10) e `s4_original_tol.jl` (`SCRATCH`), este com uma **cópia** do original alterada só no `init()`.

| Medida | Resultado | Fonte |
|---|---|---|
| Multiplicadores de Floquet dos 4 atratores do HD | \|λ\| = 0,7645 (período 1), igual a e^(−aT/2) | `[diag:tol]` |
| SCM HD 200×222: células com outro destino (padrão × Vern9 1e-8) | 62,10%; `compare_basins` = 0,8306 | `[diag:tol]` |
| SCM HD 1000×1000: idem | 81,51%; `compare_basins` = 0,7881 | `[diag:tol]` |
| Convergência Vern9 1e-8 × 1e-10 (células com outro destino) | 0,0541% (200×222); 0,1892% (1000×1000) | `[diag:tol]` |
| Original padrão × cópia do original com Vern9, HD 1000×1111 | `compare_basins` = 0,8338 | `[sp:s4]` |
| SCM (com Vern9 e agrupamento) × cópia do original com Vern9 | 0,9828 (1000×1111); 0,9735 (200×222) | `[sp:s4]` |
| SC: cópia do original com Vern9 | ainda 1 atrator e 64,87% divergentes (o gatilho é a região, §4.6) | `[sp:s4]` |
| Custo do Vern9 1e-8 | ≈ 3× (Teste 1: HD 14,9 → 41,9 s; SC 84 → 250 s) | `[log:m5]`, `[log:m5b]` |
| GCM SC com Vern9 × integrador padrão | GCM × SCM 95,71% × 74,92%; fronteira 27,23% × 60,76%; Ī estável (2,6–3,0) × crescente (3,5–5,2) | `[log:m5]`, `[log:m5b]` |

**Conclusões que valem para o relatório:**

- ⚠ **Interpretação:** a concordância de ~87% entre SCM e original **não é teto do método**. Os 13% restantes vêm da tolerância 1e-3, que o próprio original usa. Com integração precisa, os dois concordam em ~98%.
- O ruído de integração é o que torna os resultados sensíveis à máquina e ao estado do integrador. Exemplos: o atrator duplicado aparecer no notebook e não no PC; os 13,3% de diferença entre o GCM com m = 1 e o SCM no 1000×1111.
- Se, no futuro, a precisão passar a ser mais importante que a fidelidade ao original, basta trocar o `init()` nos 4 pontos listados na §5. O agrupamento funciona com os dois integradores.
