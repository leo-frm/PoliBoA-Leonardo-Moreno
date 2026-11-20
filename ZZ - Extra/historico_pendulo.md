# Pêndulo amortecido forçado: novo sistema de teste do SCM e do GCM

**Data dos testes:** 19/11/2026 (varredura) e 19/11/2026 à noite (testes oficiais)
**Escopo:** um sistema novo, acrescentado aos testes M1, M4 e M5 do SCM e do GCM, sem alterar o Helmholtz-Duffing (HD), o Spring-Column (SC), o PoliBoA original nem nenhum algoritmo.

## Fontes

Todo número abaixo indica a fonte entre colchetes.

| Marcador | Arquivo (em `ZZ - Extra/`) | De onde vem |
|---|---|---|
| `[pend_scm_m1]` | `pend_scm_m1.log` | `desenvolvimento_scm/tests_milestone1.jl`, `test_scm_pendulum(200)` |
| `[pend_gcm_m1]` | `pend_gcm_m1.log` | `desenvolvimento_gcm/tests_milestone1_gcm.jl`, `test_gcm_pendulum(200)` |
| `[pend_scm_m4]` | `pend_scm_m4.log` | `desenvolvimento_scm/tests_milestone4.jl`, `test_scm_pendulum_p(200, threads=1)` e `(…, threads=4)` |
| `[pend_gcm_m4]` | `pend_gcm_m4.log` | `desenvolvimento_gcm/tests_milestone4_gcm.jl`, `test_gcm_pendulum_p` com 1 e 4 threads |
| `[pend_scm_m5]` | `pend_scm_m5.log` | `desenvolvimento_scm/tests_milestone5.jl`, bloco `▓ PÊNDULO FORÇADO` (Testes 1–5) |
| `[pend_gcm_m5]` | `pend_gcm_m5.log` | `desenvolvimento_gcm/tests_milestone5_gcm.jl`, bloco `▓ PÊNDULO FORÇADO` (Testes 1–5) |
| `[cmp]` | `pend_<scm\|gcm>_200x200_m1.png` × `…_m4.png` | comparação byte a byte com `cmp`, feita depois da M4 |
| `[pend_og_m1]` | `pend_og_m1.log` | `pend_og_m1.jl`: PoliBoA original (`ZZ - Poli BoA og/poliboa.jl`, incluído sem alteração), só o caso B, 200×200, 1 thread, `populate_basins`. Mesma configuração do `pend_sweep_og.jl` |
| `[pend_diag_atratores]` | `pend_diag_atratores.log` | `pend_diag_atratores.jl`: atratores verdadeiros por integração de alta precisão (Vern9, tol 1e-12), sem células. Só diagnóstico |
| `[pend_diag_bacias]` | `pend_diag_bacias.log` | `pend_diag_bacias.jl`: bacia de referência de alta precisão (Vern9, tol 1e-10) no centro de cada célula, e concordância do original, do SCM e do GCM com ela. Inclui os arquivos do repositório sem alterá-los. Só diagnóstico |
| `[Bucknell]` | — | notas do curso PHYS 310 da Bucknell University, https://eg.bucknell.edu/~phys310/jupyter/chaotic_pendulum.html |
| `[Baker&Gollub]` | — | Baker & Gollub, *Chaotic Dynamics: An Introduction*, 2ª ed., Cambridge University Press, 1996 (DOI 10.1017/CBO9781139170864). Fig. 3.4, pp. 46–52; páginas do exemplar enviadas pelo aluno em 19/11 (não copiadas para o repositório, por causa dos direitos autorais) |
| `[varredura]` | `pend_sweep_og.jl`, `sweep_og_B_c0.50_F1.07_W0.667.png` | varredura de 19/11 com o original (3 casos em paralelo), 200×200. Serve como figura de referência; os números estão no `[pend_og_m1]` |
| `[19/11]` | — | diagnóstico de 19/11 (script `pend_scm_threads.jl`, já removido); números informados na tarefa. O `compare_basins` = 0,767 foi reproduzido em `[pend_diag_bacias]` |
| `(derivado)` | — | calculado a partir de números citados, com a conta indicada |

Todos os testes rodaram com `julia -t 4`, em cópias temporárias dos arquivos de teste com só as chamadas do pêndulo. As cópias foram apagadas depois.

---

## 1. Sistema

**Origem:** exemplo `servomech`, seção 5.1 de https://python-control.readthedocs.io/en/latest/nonlinear.html, descrito na página como "a simple model of a spring loaded arm driven by a motor".

**Modelo da página:**

J·θ̈ = −b·θ̇ − k·r·sin θ + τ, com τ = A·cos(ω t).

Valores padrão da página: J = 100, b = 10, k = 1, r = 1.

**Forma adimensional usada:**

θ̈ + c·θ̇ + sin θ = F·cos(Ω t̃), com:
- t̃ = ω₀·t;
- ω₀ = √(kr/J) = 0,1;
- c = b/√(Jkr), F = A/(kr), Ω = ω/ω₀.

**Parâmetros:**

| Adimensional | Na página (J = 100, k = r = 1) |
|---|---|
| c = 0,5 | b = 5 |
| F = 1,07 | A = 1,07 |
| Ω = 2/3 | ω = 0,0667 |

- **Fonte dos parâmetros:** Baker & Gollub, *Chaotic Dynamics: An Introduction* (amortecimento 1/q = 0,5, ω_D = 2/3, g = 1,07) `[Baker&Gollub]`.
  - **Conferido no livro (2ª ed., 1996), cap. 3, §3.2 "Phase diagrams and Poincaré sections", Fig. 3.4:**
    - a legenda (p. 46) dá q = 2 e, para o painel (b), g = 1,07, o caso em que aparece a duplicação de período;
    - o painel (b) está na **p. 47**;
    - a figura vai de (a) a (g), nas pp. 46–52.
  - **A seção de Poincaré do painel (b) bate com o nosso atrator R2.** Ela mostra 2 pontos perto de θ ≈ −0,3 e −0,15, com ω ≈ 2,0 e 1,95. Nosso R2 é (−0,3202; 1,9912) ↔ (−0,1467; 1,9246) `[pend_diag_atratores]`.
    - O livro mostra só a órbita de uma condição inicial. O R1, coexistente, não aparece na figura.
    - ⚠ A comparação foi feita a olho, sobre a figura impressa.
  - As notas da Bucknell usam os mesmos valores, inclusive ω_d = 2/3, e citam o mesmo painel `[Bucknell]`.
  - **ω_D = 2/3 confirmado no texto:** o início da §3.2 diz que, nessa seção, ω_D = 2/3 e q = 2 ficam fixos e só a amplitude g varia. A página desse trecho não aparece na imagem enviada.
  - O livro também observa que a forma da seção de Poincaré depende da fase em que ela é tomada (Fig. 3.5). As nossas seções são tomadas em t = 0 (mod T).
  - O período 2 foi confirmado numericamente aqui `[pend_diag_atratores]` (§3.7).
  - **Citação sugerida:** BAKER, G. L.; GOLLUB, J. P. *Chaotic dynamics: an introduction*. 2. ed. Cambridge: Cambridge University Press, 1996. Fig. 3.4(b), p. 47.
- **Por que o amortecimento foi reduzido:** com o b = 10 da página (c = 1) há só 1 atrator. Fonte: exploração por simulação direta em 19/11 (script removido; registro na conversa).
- **Função no código:** `pendulum!(dx, x, (c, F, Ω), t)`, com dx[1] = x[2] e dx[2] = −c·x[2] − sin(x[1]) + F·cos(Ω·t).

**Região e `BasinProblem`:**

| Item | Valor |
|---|---|
| x₁ (θ) | [−π, π], com `pi` exato, cíclico. O 3,14 do SC criava uma faixa morta tratada como divergente |
| x₂ (θ̇) | [−3, 3] |
| Região estendida | x₁ [−π, π], x₂ [−5, 5] |
| Grid | 200×200 |
| Período de mapeamento | T = 2π/Ω = 3π |
| Transiente | 10 períodos |
| `maximum_cycles` / `maximum_extended_cycles` / `attractor_tolerance` | 1000 / 20 / 80 |
| GCM | subgrid 10×10, m = 10, sorteio `:random`, seed = 1 |
| Integrador | o mesmo do PoliBoA original (sem Vern9) |

**Esperado:** 2 atratores, ambos de **período 2** no mapa de período T, ou seja, 2 pontos cada `[pend_diag_atratores]`. As bacias verdadeiras, nos centros das células, são 0,478 / 0,522 `[pend_diag_bacias]`. Como cada método mostra esses atratores está no §3.7.

---

## 2. O que foi alterado (só acréscimos)

| Arquivo | Acréscimos |
|---|---|
| `desenvolvimento_scm/tests_milestone1.jl` | `pendulum!`, `bench_scm_pendulum`, `test_scm_pendulum`, chamada `test_scm_pendulum(200)` |
| `desenvolvimento_scm/tests_milestone4.jl` | `pendulum!`, `bench_scm_pendulum_p`, `test_scm_pendulum_p`, chamadas com 1 e 4 threads |
| `desenvolvimento_scm/tests_milestone5.jl` | `pendulum!`, `make_pendulum_bp`, bloco `▓ PÊNDULO FORÇADO` (Testes 1–5), linha `[PEND]` no "Resumo Final" |
| `desenvolvimento_gcm/tests_milestone1_gcm.jl` | `pendulum!`, `bench_gcm_pendulum`, `test_gcm_pendulum`, chamada `test_gcm_pendulum(200)` |
| `desenvolvimento_gcm/tests_milestone4_gcm.jl` | `pendulum!`, `bench_gcm_pendulum_p`, `test_gcm_pendulum_p`, chamadas com 1 e 4 threads |
| `desenvolvimento_gcm/tests_milestone5_gcm.jl` | `pendulum!`, `make_pendulum_bp` (métodos `Int` e `Vector`), `make_pendulum_gp`, bloco `▓ PÊNDULO FORÇADO` (Testes 1–5), linhas `[PEND]` no "Resumo Final" |

- **Sem Teste 6 para o pêndulo.**
- **Inalterados:** HD, SC, `plotbasins`, `milestone2_gcm.jl` e o PoliBoA original.
- **`git diff`:** 297 inserções e 4 "remoções". As remoções são as últimas linhas antigas, que só ganharam a quebra de linha final.
- **Figuras:** salvas com `plot!(plt, dpi=300)` em `ZZ - Extra/pend_<scm|gcm>_200x200_<m1|m4>.png`.

---

## 3. Resultados

### 3.1 Atratores e bacias (grid 200×200)

| Método | Atratores | Bacias (células) | Frações (derivado: células / 40 000) | Pontos do atrator | Fonte |
|---|---|---|---|---|---|
| SCM | **2** | 19 443 / 20 557 | 0,486 / 0,514 | #1: 2 células em (−2,466; 1,155) e (−2,215; 1,335); #2: 4 células perto de (−0,14 a −0,30; 1,94 a 2,00) | `[pend_scm_m1]` |
| GCM | **2** | 19 814 / 20 186 | 0,495 / 0,505 | #1: 9 células, de (−2,466; 1,155) a (−2,246; 1,305); #2: 7 células, de (−0,298; 1,995) a (−0,141; 1,935) | `[pend_gcm_m1]` |
| Original | **2** | 20 486 / 19 514 | 0,512 / 0,488 | #1: 1 célula em (−2,466; 1,155); #2: 1 célula em (−0,141; 1,935). O original só lista a célula mais visitada (§3.7) | `[pend_og_m1]` |
| Referência (Vern9) | **2**, período 2 | 19 119 / 20 881 | 0,478 / 0,522 | R1: (−2,2077; 1,3461) ↔ (−2,4770; 1,1494); R2: (−0,3202; 1,9912) ↔ (−0,1467; 1,9246) | `[pend_diag_atratores]`, `[pend_diag_bacias]` |

**Análise do GCM** `[pend_gcm_m5]` (Teste 2):
- 3 grupos persistentes, incluindo o sumidouro, com períodos 1, 1 e 1. O período 1 vem das transições entre células vizinhas dentro do grupo; a dinâmica verdadeira tem período 2 (§3.7);
- 16 células persistentes (0,04%);
- domicílio único: 14 904 células (37,26%);
- domicílio múltiplo: 25 096 células (62,74%).

### 3.2 M4: 1 thread × 4 threads

| Método | Imagem da M1 × imagem da M4 (4 threads) | Fonte |
|---|---|---|
| GCM | **idênticas byte a byte** | `[cmp]` |
| SCM | **diferentes** (primeiro byte diferente: 89) | `[cmp]` |

- **Por que o SCM difere:** o SCM não tem `reset_dt`, então o resultado depende da ordem em que as células são integradas.
- **Tamanho da diferença** `[19/11]`: com 4 threads, 26 346 de 40 000 destinos mudaram em relação ao serial; `compare_basins` = 0,767; ainda 2 atratores.
- **Efeito na exatidão** `[pend_diag_bacias]`: nenhum. Contra a referência, o SCM acerta 72,02% com 1 thread e 71,96% com 4. Pelo menos 94,5% das diferenças entre as duas execuções estão em células de fronteira (derivado, §3.7).
- Os logs da M4 não imprimem o número de atratores.

### 3.3 M5, Teste 2: concordância de bacias

| Comparação | Grid inteiro | Domicílio único (14 904 células) | Domicílio múltiplo (25 096 células) | Fonte |
|---|---|---|---|---|
| SCM × Original | 71,37% | — | — | `[pend_scm_m5]` |
| GCM × Original | 70,63% | 97,01% | 54,96% | `[pend_gcm_m5]` |
| GCM × SCM | 70,36% | 96,25% | 54,98% | `[pend_gcm_m5]` |

### 3.4 M5, Testes 1 e 3: tempo e custo (1 thread, 3 repetições)

| Medida | SCM | GCM | Original | Fonte |
|---|---|---|---|---|
| Teste 1: tempo médio | 1,7198 ± 0,0277 s | 17,9875 ± 0,2542 s | — | `[pend_scm_m5]`, `[pend_gcm_m5]` |
| Teste 1: bytes alocados | 28,05 MB | 506,82 MB | — | idem |
| Teste 1: integrações | 40 000 | 400 000 | — | idem |
| Teste 3: tempo médio | 1,7254 s | 17,6642 s | 1,5008 s (log do SCM) / 1,5152 s (log do GCM) | idem |

**Razões do Teste 3:**

| Comparação | Tempo | Bytes alocados | Integrações | Fonte |
|---|---|---|---|---|
| SCM / Original | 1,15× | 1,28× | 1,00× | `[pend_scm_m5]` |
| GCM / Original | 11,66× | 23,04× | 10,00× | `[pend_gcm_m5]` |
| GCM / SCM | 10,02× | 18,07× | 10,00× | `[pend_gcm_m5]` |

### 3.5 M5, Teste 4: memória

| Grid | SCM: MB / MB/kCell | GCM: MB / MB/kCell | GCM: Ī | GCM: mapa (MB) | Fonte |
|---|---|---|---|---|---|
| 100×100 | 6,7 / 0,667 | 114,9 / 11,494 | 6,311 | 0,521 | `[pend_scm_m5]`, `[pend_gcm_m5]` |
| 200×200 | 26,7 / 0,669 | 483,3 / 12,084 | 7,134 | 2,335 | idem |
| 400×400 | 107,1 / 0,669 | 1942,0 / 12,138 | 7,838 | 10,197 | idem |

⚠ **Interpretação:** o MB/kCell do SCM é constante, então a memória é O(N). No GCM ele sobe pouco (11,5 → 12,1), acompanhando o Ī, que cresce com N (6,3 → 7,8). Ou seja, o número médio de imagens por célula aumenta com a resolução. É o mesmo comportamento do SC com o integrador do original (`historico_GCM.md`, §2.4).

### 3.6 M5, Teste 5: speedup (grid 200×200)

| Threads | SCM | GCM | Fonte |
|---|---|---|---|
| 1 | 1,7304 s (1,00×) | 17,9986 s (1,00×) | `[pend_scm_m5]`, `[pend_gcm_m5]` |
| 2 | 0,8591 s (2,01×) | 8,9897 s (2,00×) | idem |
| 4 | 0,4773 s (**3,63×**) | 5,5615 s (**3,24×**) | idem |

### 3.7 Diagnóstico: referência de alta precisão (fora dos testes oficiais)

Feito para resolver as pendências: os ~71% do Teste 2, o SCM com M1 ≠ M4 e o período citado do Baker & Gollub.
- Vern9 só entra como referência; os métodos continuam com o integrador do original.
- Nada do repositório foi alterado.
- Grid 200×200, `julia -t 4`.

**Atratores verdadeiros** `[pend_diag_atratores]`: 2, ambos de período 2.
- R1: (−2,2077; 1,3461) ↔ (−2,4770; 1,1494);
- R2: (−0,3202; 1,9912) ↔ (−0,1467; 1,9246);
- de 1000 condições iniciais aleatórias, 465 vão para R1 e 535 para R2. Não há outro atrator.

**Como cada método enxerga esses atratores:**
- **Integrador do original (reltol 1e-3):** a órbita alterna entre os dois pontos de R1, mas com desvios de até ~0,07 em θ (≈ 2 células) `[pend_diag_bacias]`.
- **SCM:** o atrator #1 tem exatamente as 2 células dos dois pontos de R1 `[pend_scm_m1]`.
- **GCM:** os grupos persistentes cobrem os dois pontos (9 e 7 células) `[pend_gcm_m1]`. O "período 1" do grupo vem das transições entre células vizinhas, não da dinâmica.
- **Original:** só lista as células que chegaram a `attractor_tolerance` visitas antes de a busca parar. Com os desvios do integrador, sobra 1 célula por atrator `[pend_og_m1]`.

**Bacia de referência** `[pend_diag_bacias]`:
- o centro de cada célula é integrado até ficar a menos de 1e-3 de um ponto de R1 ou R2 (mediana de 19 períodos, máximo 117; todas convergiram);
- resultado: R1 com 19 119 células (0,478) e R2 com 20 881 (0,522);
- **29 307 células (73,27%) são de fronteira**: na vizinhança 3×3 (x₁ cíclico), têm um vizinho com outro rótulo. As outras 10 693 são interiores.

| Método | Atratores | Bacias | × referência: total | Interior | Fronteira |
|---|---|---|---|---|---|
| Original, 1 thread | 2 | 0,512 / 0,488 | 71,98% | 98,21% | 62,40% |
| SCM, 1 thread | 2 | 0,486 / 0,514 | 72,02% | 97,52% | 62,72% |
| SCM, 4 threads | 2 | 0,479 / 0,521 | 71,96% | 97,73% | 62,56% |
| GCM, 4 threads | 2 | 0,495 / 0,505 | 71,79% | **99,08%** | 61,83% |
| Original, 4 threads | **4** | 0,125 / 0,494 / 0,247 / 0,134 | 53,10% | 69,03% | 47,29% |
| SCM + `reset_dt` (variante só do diagnóstico, **não aplicada**), 1 ou 4 threads | **3** | 0,485 / 0,312 / 0,203 | 58,33% | 75,40% | 52,10% |

A concordância entre execuções, medida com `compare_basins`, é esta `[pend_diag_bacias]`:

| Comparação | `compare_basins` |
|---|---|
| Original × SCM | 0,7137 |
| Original × GCM | 0,7063 |
| SCM × GCM | 0,7036 |
| SCM, 1 × 4 threads | 0,7671 |
| SCM + `reset_dt`, 1 × 4 threads | 1,0000 |
| Original, 1 × 4 threads | 0,5530 |

As três primeiras linhas reproduzem os valores do Teste 2 (§3.3).

**Conclusões:**
1. **Os ~71% não são defeito de nenhum método.**
   - Os três métodos ficam igualmente perto da referência (71,8% a 72,0%).
   - No interior das bacias, acertam de 97,5% a 99,1%.
   - Na fronteira, que cobre 73% do grid a 200×200, todos ficam perto de 62%, pouco acima de cara ou coroa.
   - A concordância baixa entre eles vem só da fronteira.
2. **O SCM com M1 ≠ M4 não muda a exatidão.**
   - As células interiores são 10 693. Contra a referência, o SCM erra no máximo 265 delas com 1 thread (2,48%) e 243 com 4 threads (2,27%). Logo, as duas execuções divergem em no máximo 508 células interiores.
   - No total, as duas execuções divergem em 40 000 × (1 − 0,7671) ≈ 9 316 células. Assim, pelo menos (9 316 − 508) / 9 316 = 94,5% das diferenças estão na fronteira (derivado).
3. **Pôr `reset_dt` no SCM não é uma correção simples.** A variante fica reprodutível (1,0000), mas neste sistema cria um terceiro atrator e cai para 58%. Não foi aplicada: mudaria o algoritmo do SCM e as figuras do HD e do SC. A decisão é do aluno.
4. **O original também depende das threads.** Com 4 threads, encontra 4 atratores. Os testes M5 rodam o original com 1 thread, então os números do §3.3 não são afetados.

---

## 4. Resumo e observações

- **Número de atratores:** 2 no SCM, no GCM e no original, como esperado. Os dois são de período 2, confirmado com integração de alta precisão (§3.7).
- **Concordância:**
  - os três métodos concordam em ~70–71% no grid inteiro;
  - nas células de domicílio único, o GCM concorda 97,0% com o original e 96,3% com o SCM;
  - as discordâncias ficam na fronteira entre bacias. No GCM, as células de domicílio múltiplo cobrem 62,74% do grid; na referência, 73,27% das células são de fronteira;
  - **Verificado** `[pend_diag_bacias]`: contra uma referência de alta precisão, os três métodos acertam ~72% no total e 97,5% a 99,1% no interior (§3.7). A ver também `sweep_og_B_c0.50_F1.07_W0.667.png`.
- **Custo:** GCM ≈ 10× o SCM em tempo e em integrações (m = 10). O SCM leva ≈ 1,15× o tempo do original.
- **Paralelização:**
  - o GCM é reprodutível entre 1 e 4 threads (M1 = M4);
  - o SCM não é, pela falta de `reset_dt` (§3.2), mas a exatidão é a mesma com 1 e 4 threads (§3.7);
  - pôr `reset_dt` no SCM piorou o resultado neste sistema (3 atratores). Não foi aplicado; nenhum algoritmo foi alterado.
- **Teste 6:** não foi feito para o pêndulo, por decisão do aluno.
