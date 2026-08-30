# Revisão técnica — ElectionForensics.jl

**Revisor:** auditoria de engenharia de pacotes Julia + estatística forense eleitoral
**Commit:** `590903f` (+ `test/runtests.jl` modificado no working tree)
**Data:** 2026-08-29 · **Julia:** 1.12.7 (também verificado em 1.10 LTS e 1.6.7)
**Scripts de evidência:** `validation/` · **Patches:** `validation/correcoes.patch`

---

## 1. Veredito

> ## `NÃO APTO` para release

Dois dos quatro testes da bateria são **anti-conservadores por construção**: no
DGP eleitoral mais comum no Brasil (seções de tamanho homogêneo), o teste do
penúltimo dígito rejeita eleições limpas em **50–100 %** dos casos e o teste de
Benford 2BL em **72 %** — e ambos são executados por default em
`forensics_report`. Um pacote de forense eleitoral que acusa fraude em dados
limpos com essa frequência não pode ser publicado no General Registry.

A boa notícia: **a engenharia está correta e as correções são localizadas.** A
maquinaria χ², as probabilidades de Benford, os limiares de Nigrini, a fórmula
do p-valor Monte Carlo e o teste de Rozenas conferem com a literatura e com
implementações independentes. Os patches em `validation/correcoes.patch` já
foram aplicados e verificados: 105/105 testes passam, `Aqua.test_all` limpo, e a
taxa de erro tipo I da 2BL cai de 0,72 para 0,058.

### Escopo: o briefing não bate com o repositório

O pedido de revisão menciona modelo de Klimek et al. (fingerprint 2D),
quadratura de Gauss–Hermite, resampling bayesiano de Rozenas, bootstrap
BCa, paralelização por threads, extensões Makie, tema Dracula e retornos
Tables.jl. **Nada disso existe no código.** O pacote tem 521 linhas em 5 arquivos
e implementa Benford (1BL/2BL), último/penúltimo dígito e um teste de frações
coarse. Os itens correspondentes do briefing estão marcados como
*não aplicável* na seção 7, não como aprovados.

---

## 2. Inventário: método → arquivo:linha → referência

| Método | Local | Referência alegada | Confere? |
|---|---|---|---|
| Probabilidades Benford 1º dígito | `src/benford.jl:11` | Benford (1938) | ✅ exato |
| Probabilidades Benford 2º dígito | `src/benford.jl:13` | Mebane (2008) | ✅ exato |
| Limiares MAD de conformidade | `src/benford.jl:20-23` | Nigrini (2012), Tab. 5.1 | ✅ exato |
| χ² de aderência + p-valor | `src/benford.jl:85-87` | — | ✅ exato, mas null inadequado |
| Uniformidade do último dígito | `src/lastdigit.jl:51-55` | Beber & Scacco (2012) | ✅ calibrado |
| Uniformidade do penúltimo dígito | `src/lastdigit.jl:38` | Beber & Scacco (2012) | ❌ null falso |
| Frações de Farey (coarse) | `src/rozenas.jl:13` | Rozenas (2017) | ✅ |
| Nula por bootstrap paramétrico | `src/rozenas.jl:99-112` | Rozenas (2017), §3 | ⚠️ aproximação, não o RKD do artigo |
| p-valor Monte Carlo | `src/rozenas.jl:119` | Davison & Hinkley (1997) | ✅ convenção correta |
| Banda de Silverman | `src/utils.jl:67-76` | Silverman (1986), eq. 3.31 | ✅ |
| Bateria integrada | `src/report.jl:28-63` | — | ❌ inclui teste inválido |

Exportações: 12 símbolos (4 funções de dígito, 3 testes, `coarse_fractions`,
`forensics_report`, 3 tipos de resultado). Sem extensões, sem `weakdeps`,
sem `docs/`, sem CI.

---

## 3. O que está CORRETO (verificado por execução)

Não é um relatório punitivo. Estes pontos foram checados e passam:

| Item | Evidência |
|---|---|
| **Probabilidades de Benford** | 1BL `[.30103, .176091, …, .045757]` e 2BL `[.119679, …, .084997]` — batem com a literatura em todas as casas; somam 1,0 |
| **Limiares de Nigrini** | `(0.006, 0.012, 0.015)` para 1BL e `(0.008, 0.010, 0.012)` para 2BL — idênticos à Tab. 5.1 de Nigrini (2012) |
| **Graus de liberdade** | 8 (1BL) e 9 (2BL, último dígito) — corretos, nenhum parâmetro estimado |
| **Maquinaria χ² exatamente calibrada** | Com dígitos sorteados da própria lei de Benford: KS p = 0,29–0,91; rejeição a α=0,05 → **0,046–0,054** (n=500 e 2000, erro padrão 0,003). `validation/01b` |
| **Referência cruzada com scipy** | 12 comparações de χ², p-valor e MAD sobre datasets idênticos: discrepância relativa máxima **2,8 × 10⁻¹²**. `validation/10b` |
| **Último dígito bem calibrado** | Erro tipo I = 0,044–0,058 para seções de 60 a 5 000 eleitores, e estável em n de 500 a 30 000. `validation/02` |
| **Poder do último dígito** | Arredondamento decimal em 5 % das seções → poder **1,000**. `validation/06` |
| **Rozenas: erro tipo I ≤ α** | 12 cenários (Beta(8,6), Beta(2,2), Beta(1.2,8); totais 20–120 a 150–900): rejeição 0,016–0,055 a α=0,05. Nunca inflado. `validation/03` |
| **Rozenas: poder real** | Metas redondas em 2 % das seções → 0,47; em 5 % → **0,99**. `validation/06` |
| **A nula de Rozenas não distorce a variância** | sd da nula / sd verdadeiro de T = 0,87–0,99 em 4 regimes. `validation/04` |
| **`h` é robusto** | Variar a banda em 32× (h/8 a 4h) move o p-valor de 0,088 a 0,135. `validation/04` |
| **Aritmética racional exata** | `votes[i] // totals[i]` compara sem erro de float — é o primitivo correto para o método |
| **Fórmula do p-valor MC** | `(1 + #{T* ≥ T}) / (B + 1)` — convenção correta, inclui empates (conservadora) |
| **Reprodutibilidade** | Mesmo `rng` ⇒ mesmo p-valor. Pacote é single-thread: **não há dependência do número de threads** (o risco levantado no briefing não se aplica) |
| **Aqua** | Sem ambiguidades, sem pirataria de tipo, sem deps obsoletas, sem exports indefinidos. Só falha em `[compat]` |
| **JET** | `report_package` não aponta nada em código do pacote (só internals de `Base`/`Statistics`) |
| **Compat `julia = "1.6"`** | Verificado: 76/76 testes passam em Julia 1.6.7 e em 1.10 LTS |
| **Tempo de carga** | 0,38 s; TTFX 0,77 s (Benford), 1,47 s (Rozenas). Sem problema |
| **Precompilação** | Sem warnings |

Também testei uma correção alternativa que **não** foi necessária: substituir o
null uniforme do **último** dígito por um null reamostrado não muda nada
(0,0675 → 0,0625 sob H0; poder idêntico). O teste do último dígito está bom
como está.

---

## 4. Tabela de achados

### 🔴 Crítico

| # | Local | Problema | Evidência | Correção |
|---|---|---|---|---|
| **C1** | `src/lastdigit.jl:33-60`, `src/report.jl:53` | O teste do **penúltimo dígito** usa null Uniforme{0..9}, que é falso para contagens eleitorais. Roda por default em `forensics_report`. | Eleições **limpas** simuladas (m=2000): rejeição a α=0,05 de **0,503** (seções 150–900), **1,000** (100–200), **0,667** (Beta(1,1)). Num dataset limpo tipo TSE o `forensics_report` atual devolve **p = 4,8 × 10⁻⁸**. Sob H₀ a distribuição do penúltimo dígito é decrescente (10,95 % → 8,87 %), não uniforme. `validation/07` | Tornar opt-in (`penultimate = false` no relatório), campo `null_valid`, aviso e docstring. Patch aplicado. |
| **C2** | `src/benford.jl:87`, `src/report.jl:45` | O null da **lei de Benford** não descreve contagens eleitorais limpas; o p-valor não controla erro tipo I. Default `digit = 2` é apresentado como a escolha segura. | Rejeição em eleições limpas: 1BL = **1,000** em 6 de 8 DGPs; 2BL = **0,710** (seções 200–400, caso TSE), 0,332 (150–900), 0,218 (Pareto). Só fica em α quando os tamanhos de seção são log-normais com σ grande. `validation/01, 01b, 12` | `null = :resampled` como default (jitter gaussiano em `log10`), mantendo `:benford` como descritivo. **Verificado: 0,710 → 0,058.** Patch aplicado. |

### 🟠 Alto

| # | Local | Problema | Evidência | Correção |
|---|---|---|---|---|
| **A1** | `Manifest.toml` (versionado) | Manifest de pacote no repositório. Foi resolvido com Julia 1.12.6 e **quebra `Pkg.instantiate()` na versão mínima declarada**. | Julia 1.6.7 + este Manifest: `ERROR: AssertionError: sourcepath !== nothing`. Sem o Manifest, 76/76 passam. | `git rm --cached Manifest.toml` + `.gitignore`. Patch aplicado. |
| **A2** | `Project.toml:12-14, 16-17` | Faltam `[compat]` para `Printf`, `Random`, `Statistics` e `Test`. Bloqueia AutoMerge do General. | `Aqua.test_deps_compat` falha em `deps` e `extras`. | Adicionar `= "1"` a cada; `Test = "1"`. Patch aplicado. |
| **A3** | `README.md:33` | Exemplo do README **não executa**: `max_denom` não é keyword de `rozenas_test`. | `MethodError: no method matching rozenas_test(::Vector{Int64}, ::Vector{Int64}; max_denom, B)` | Trocar por `fractions = coarse_fractions(10)`. Patch aplicado. |
| **A4** | `src/rozenas.jl:132-143` | O `show` imprime z-scores por fração **sem correção de multiplicidade**. São máximos de 31 estatísticas: o usuário lê ruído como evidência. | Sob H₀: maior z tem mediana **2,41**, q90 = 3,85; `P(algum z > 2) = 0,678`; em média **1,14** frações com z > 2 por eleição limpa. `validation/11` | Campo `qvalues` (Benjamini–Hochberg) e exibi-lo na tabela. **Verificado: 0 falsas descobertas em 200 eleições limpas.** Patch aplicado. |
| **A5** | `src/rozenas.jl:102` | `clamp(share + h·randn, 0, 1)` empilha massa **exatamente** em 0 e 1, onde a Binomial degenera e nunca produz fração coarse — artefato não controlado na nula. | Fração de sorteios colados na fronteira: 0,7 % (Beta(2,2)), 13,9 % (Beta(.5,.5)), **23,8 %** (Beta(.3,.3)). Deflaciona E[T*] de 18,2 para 15,9. `validation/05` | Reflexão em vez de truncamento. Patch aplicado. |
| **A6** | repositório | Sem `.github/workflows`, sem `docs/`, sem TagBot, sem CompatHelper, sem `.gitignore`. | `ls .github docs .gitignore` → inexistentes | CI (1.10/1/pre × 3 SOs), TagBot, CompatHelper, doctests. Patch aplicado. |
| **A7** | `src/benford.jl:33`, `lastdigit.jl:3`, `rozenas.jl:17` | Os 3 tipos exportados **não têm docstring**, mas são referenciados com `[`X`](@ref)`. Documenter com `checkdocs = :exports` falharia e os links quebram. | Verificado por introspecção: `BenfordResult, LastDigitResult, RozenasResult` sem documentação. | Docstrings adicionadas. Patch aplicado. |

### 🟡 Médio

| # | Local | Problema | Evidência |
|---|---|---|---|
| M1 | `src/benford.jl:85`, `lastdigit.jl:53` | Nenhuma checagem da **regra de Cochran** (casela esperada ≥ 5). O χ² é reportado como válido com n=1. | `benford_test([10]; digit=2)` devolve `p = 0.600` com um único dado. Corrigido no patch (aviso `min_expected`). |
| M2 | `test/runtests.jl` (todo) | Asserções são limiares sobre **um sorteio**, com um `MersenneTwister` compartilhado sequencialmente: inserir um teste no meio muda todos os seguintes. | `@test r.pvalue > 0.05` com B=199 sobre dados cujo p verdadeiro é ~0,11: sd do p-valor MC = **0,028**. Ver §6. |
| M3 | `src/lastdigit.jl:35` | `min_value = 10` é permissivo demais para `:last`. | Seções de 30 eleitores: erro tipo I = **0,297**. Mistura realista (N ∈ 10..400, n=10 000): 0,122. `validation/02` |
| M4 | assinaturas | `AbstractVector{<:Integer}` rejeita `Vector{Float64}` e `Vector{Union{Missing,Int}}` com `MethodError` cru — formatos comuns de CSV eleitoral. | `validation/09` |
| M5 | `src/rozenas.jl:99-112` | A nula adiciona ruído binomial **sobre shares que já contêm ruído binomial** (Var = σ²+h²+ruído, duas vezes). Torna o teste conservador onde os shares se acumulam nas bordas. | Beta(.3,.3): E[T*] = 18,2 vs E[T_obs] = 13,8; rejeição **0,0025** contra α=0,05. Remédio principiado: bootstrap suavizado com correção de variância (Silverman & Young 1987). `validation/05` |
| M6 | `src/rozenas.jl:63` | Sem guarda de `m` mínimo: `rozenas_test([50],[100])` roda e devolve p = 0,08. | `validation/09`. Aviso adicionado no patch. |
| M7 | `src/utils.jl:75` | `_silverman` cai num literal mágico `0.01` quando a dispersão é nula, sem sinalizar. | `rozenas_test(fill(50,500), fill(100,500))` → `h = 0.01` silenciosamente. |
| M8 | `README.md:39`, `src/rozenas.jl:40` | O método é rotulado "teste de Rozenas (2017)". É uma **aproximação frequentista**; o artigo usa um estimador de densidade kernel reamostrada bayesiano. A docstring é honesta ("aproximação"), o README não qualifica. | Comparação com Rozenas (2017), §3. |
| M9 | tipos de retorno | Nenhum resultado é compatível com Tables.jl, contrariando a convenção do projeto. Esboço: `Tables.istable(::Type{BenfordResult}) = true` + `rowtable` sobre `(digit, count, observed, expected)`. | — |
| M10 | `src/report.jl:53` | `forensics_report` lança `ArgumentError` se todas as contagens forem < 100 (eleições municipais pequenas). | Corrigido indiretamente: penúltimo saiu do default. |
| M11 | `src/benford.jl:88` + `:show` | MAD e χ² podem discordar (MAD ignora n) e o `show` exibe ambos sem explicar. | Dados arredondados: χ² p = 2 × 10⁻⁴ mas conformidade Nigrini = `:marginal`. |

### 🟢 Baixo

| # | Local | Problema |
|---|---|---|
| B1 | `src/benford.jl:27-30` | `_nigrini_conformity` usa `<` estrito; os intervalos de Nigrini são fechados no limite superior. |
| B2 | `src/benford.jl:75` | `getd = digit == 1 ? first_digit : second_digit` cria uma `Union` de tipos de função (instabilidade). |
| B3 | `resultados_forensics/*.png` | 3 PNGs versionados, não referenciados pelo README nem pelo código. |
| B4 | `src/benford.jl:72`, `lastdigit.jl:44` | `@warn` não suprimível; polui laços de simulação. (Resolvido pelo `warn` keyword do patch.) |
| B5 | `src/rozenas.jl:105` | `y // totals[i]` num laço quente: `hash(::Rational)` passa por `ldexp`/`reinterpret` (apontado pelo JET). Trocar por chave `(num, den)` normalizada evitaria o custo. |
| B6 | `src/utils.jl:34` | `last_digit_test` com `min_value` negativo aceita contagens negativas e usa `abs` silenciosamente. |

---

## 5. Patches

Todos os patches de severidade Crítica e Alta estão em
**`validation/correcoes.patch`** (1 010 linhas, 13 arquivos). Foram aplicados a
uma cópia limpa e verificados:

```
Test Summary:        | Pass  Total   Time
ElectionForensics.jl |  105    105  25.4s
     Testing ElectionForensics tests passed
```

`Aqua.test_all(ElectionForensics)` passa **sem exclusões** (era 2 falhas).

### Verificação empírica dos patches

| Métrica | Antes | Depois |
|---|---|---|
| 2BL, erro tipo I, seções 200–400 (TSE) | **0,710** | **0,058** |
| 2BL, erro tipo I, seções 150–900 | 0,332 | 0,038 |
| 2BL, erro tipo I, totais LN(6, 1.2) | 0,062 | 0,048 |
| Penúltimo dígito no `forensics_report` | roda por default | opt-in + `null_valid` |
| Rozenas, poder (metas redondas, ε=5 %) | 0,988 | **0,988** (preservado) |
| Rozenas, poder (ε=2 %) | 0,472 | 0,510 |
| Frações com q < 0,05 em eleições limpas | 1,14 com z > 2 | **0,000** |

### Resumo dos diffs

**`src/benford.jl`** — introduz `Base.@kwdef struct BenfordConfig` (`digit`,
`null`, `B`, `log_bandwidth`, `min_n`, `min_expected`, `warn`), o null
reamostrado `_resampled_pvalue`, o campo `pvalue_asymptotic` e a checagem de
Cochran:

```diff
+Base.@kwdef struct BenfordConfig
+    digit::Int = 2
+    null::Symbol = :resampled
+    B::Int = 999
+    log_bandwidth::Float64 = 0.05
+    min_n::Int = 100
+    min_expected::Float64 = 5.0
+    warn::Bool = true
+end
+
+# null reamostrado: jitter gaussiano em log10 apaga a estrutura do dígito
+# preservando o formato macro da distribuição de contagens.
+function _resampled_pvalue(v, digit, expected, chi_obs, cfg, rng)
+    ...
+            y = round(Int, 10^(lv[rand(rng, 1:n)] + cfg.log_bandwidth * randn(rng)))
```

**`src/lastdigit.jl`** — `Base.@kwdef struct DigitTestConfig`,
`_orders_of_magnitude`, campo `null_valid` e aviso:

```diff
+_orders_of_magnitude(v) =
+    log10(quantile(Float64.(v), 0.95) / max(quantile(Float64.(v), 0.05), 1.0))
...
+    null_valid = cfg.position === :last ? true : oom ≥ 2.0
+        if cfg.position === :penultimate && !null_valid
+            @warn "null uniforme do penúltimo dígito NÃO é válido para estes \
+                dados: contagens cobrem $(round(oom, digits = 2)) décadas (< 2)."
```

**`src/rozenas.jl`** — reflexão na fronteira e q-valores BH:

```diff
-            p = clamp(shares[i] + hval * randn(rng), 0.0, 1.0)
+            # reflexão nas fronteiras: `clamp` empilharia massa exatamente em
+            # 0 e 1, onde Binomial degenera e nunca produz fração coarse — com
+            # shares em U (eleições polarizadas) isso atingia 22 % dos sorteios
+            # e distorcia a nula (validation/05_rozenas_clamp.jl).
+            p = shares[i] + hval * randn(rng)
+            while p < 0.0 || p > 1.0
+                p < 0.0 && (p = -p)
+                p > 1.0 && (p = 2.0 - p)
+            end
...
+    praw = [(1 + count(≥(observed[j]), @view null_counts[:, j])) / (B + 1) for j in 1:J]
+    qvalues = _bh_adjust(praw)
```

**`src/report.jl`** — `penultimate::Bool = false`, `benford_null = :resampled`,
validação de entrada e um bloco `!!! warning` sobre o que a bateria **não**
detecta.

**`src/utils.jl`** — `_bh_adjust` (Benjamini–Hochberg, com monotonicidade
forçada e testado contra valores analíticos).

Aplicar com:

```bash
git rm --cached Manifest.toml
git apply validation/correcoes.patch
julia --project=. -e 'using Pkg; Pkg.test()'
```

---

## 6. Eixo D — a suíte de testes

A suíte original (76 asserções) tem boa cobertura de **caminhos** e de
**validação de argumentos**, mas quase nenhuma de **valores esperados
analíticos**. Ela verifica sobretudo que as funções não quebram e que
p-valores caem do lado certo de um limiar.

**O que já existe e é bom:** `@test_throws` para todos os erros de argumento;
verificação de que `sum(expected) ≈ 1`; testes de `n_excluded`; teste de
fraude injetada com poder alto.

**O que falta:**

1. **Valores analíticos.** Nenhum teste compara `benford_expected(1)[1]` com
   `log10(2) = 0.30103`, nem `chi2` com um valor calculado à mão. Um dataset
   fixo de 10 elementos com χ² conferido em papel valeria mais que os 5 testes
   de "não quebra".
2. **Afirmações sobre taxa, não sobre um sorteio.** `@test r.pvalue > 0.05` com
   B=199 sobre dados cujo p verdadeiro é ~0,11 cruza o limiar por acaso: o sd
   Monte Carlo do p-valor é **0,028** (medido em `validation/04`). O patch
   substitui isso por rejeições contadas sobre 40 réplicas com seeds fixas —
   escrevi o teste da forma errada primeiro e ele falhou, o que é a
   demonstração do ponto.
3. **RNG encadeado.** Um único `MersenneTwister(2026)` é consumido
   sequencialmente por todos os testsets. Inserir um teste no meio altera os
   dados de todos os posteriores. Use uma semente por testset.
4. **Nenhum teste de calibração.** Nenhuma asserção sobre erro tipo I. Os dois
   defeitos críticos deste relatório teriam sido pegos por um único teste de
   1 000 réplicas sob H₀.
5. **Sem Aqua**, sem doctests, sem cobertura medida.
6. **Casos-limite ausentes:** n=1, todas as contagens iguais, `totals = 1`,
   contagens gigantes, `missing`.

O patch leva a suíte a 105 asserções, adicionando: valores analíticos para
`_bh_adjust`, o teste de calibração por taxa da 2BL, a regressão do
`null_valid`, e `Aqua.test_all`.

---

## 7. Checklist de release

**Project.toml e metadados**
- [x] `name`, `uuid`, `version` presentes e coerentes
- [x] `[compat]` para `Distributions` e `julia`
- [ ] → [x] `[compat]` para `Printf`, `Random`, `Statistics`, `Test` *(patch)*
- [x] `[extras]` / `[targets]` de teste declarados
- [n/a] `[weakdeps]` / `[extensions]` — o pacote não tem extensões
- [x] `julia = "1.6"` **verificado** (76/76 em 1.6.7) — o patch sobe para `"1.10"` (LTS atual)
- [ ] → [x] `Manifest.toml` fora do controle de versão *(patch)*
- [ ] → [x] `.gitignore` *(patch)*

**Documentação**
- [x] `LICENSE` presente (MIT)
- [x] README com seção de métodos e referências
- [ ] → [x] README com exemplo **executável** *(patch: `max_denom` corrigido)*
- [x] Docstrings nas 9 funções exportadas
- [ ] → [x] Docstrings nos 3 tipos exportados *(patch)*
- [ ] `docs/` com Documenter — **ausente**, não corrigido pelo patch
- [ ] Doctests — nenhum docstring tem bloco `jldoctest`
- [ ] CHANGELOG

**CI**
- [ ] → [x] Matriz Julia (1.10 LTS, 1, pre) *(patch)*
- [ ] → [x] Múltiplos SOs (Linux, macOS, Windows) *(patch)*
- [ ] → [x] TagBot, CompatHelper *(patch)*
- [ ] → [x] Cobertura (Codecov) *(patch)*
- [ ] Job de doctests configurado mas sem `docs/` para rodar

**Qualidade**
- [x] `Aqua`: ambiguidades, pirataria, deps obsoletas, exports — limpo
- [ ] → [x] `Aqua`: `[compat]` *(patch)*
- [x] `JET`: nada em código do pacote
- [x] `] test` passa do zero em ambiente limpo (1.6.7, 1.10.12, 1.12.7)
- [x] Precompila sem warnings
- [x] Tempo de carga aceitável (0,38 s)
- [ ] → [x] `@warn` suprimível *(patch: keyword `warn`)*

**Correção estatística**
- [x] Probabilidades de Benford conferidas contra a literatura
- [x] Graus de liberdade corretos
- [x] Referência cruzada com `scipy` (Δ ≤ 2,8 × 10⁻¹²)
- [x] p-valor Monte Carlo com a convenção correta
- [x] RNG reprodutível, sem dependência de threads
- [ ] → [x] Calibração sob H₀ da 2BL *(patch: 0,710 → 0,058)*
- [ ] → [x] Penúltimo dígito fora do default *(patch)*
- [ ] → [x] Multiplicidade nas frações *(patch: BH)*
- [ ] Conservadorismo residual do Rozenas sob shares em U — **não corrigido**, ver M5
- [ ] Retornos compatíveis com Tables.jl — **não implementado**, ver M9

**Não verificado**
- [ ] Dados reais do TSE — o CDN respondeu **HTTP 403** a clientes não-browser
      neste ambiente (`cdn.tse.jus.br` e `dadosabertos.tse.jus.br`). Script
      pronto em `validation/14_dados_reais_tse.jl`; rode localmente.
- [ ] Build do Documenter — não há `docs/` para compilar.

**Não aplicável (não existe no pacote)**
- Modelo de Klimek et al. / fingerprint 2D · quadratura de Gauss–Hermite ·
  resampling bayesiano de Rozenas · `logsumexp` · bootstrap percentil/BCa ·
  paralelização por threads · extensões Makie · tema Dracula

---

## 8. Versão sugerida e changelog

O pacote **nunca foi publicado** (uma única tag, `590903f`, e nenhum release no
General). Portanto **não há bump**: publique `0.1.0` já com as correções. Lançar
0.1.0 com os defeitos e corrigir em 0.2.0 colocaria uma versão que acusa fraude
em dados limpos no registro permanente.

Se preferir marcar a linha do tempo interna, use **0.1.0 → 0.2.0 (minor)**: as
mudanças são quebras de contrato (default de `null`, `forensics_report` deixa de
retornar `penultimate::LastDigitResult`, `RozenasResult` ganha um campo), mas em
0.x o minor já sinaliza incompatibilidade sob SemVer.

### CHANGELOG.md (rascunho)

```markdown
# Changelog

## [0.1.0] — não lançado

Primeira versão pública. Testes forenses de dígitos e de frações coarse.

### Adicionado
- `benford_test` (1BL/2BL) com null reamostrado (`BenfordConfig`)
- `last_digit_test` (último e penúltimo dígito, `DigitTestConfig`)
- `rozenas_test` — frações coarse com nula por bootstrap paramétrico
- `forensics_report` — bateria integrada com saída colorida
- `qvalues` (Benjamini–Hochberg) por fração em `RozenasResult`
- `null_valid` em `LastDigitResult`
- Suíte `validation/` com calibração sob H₀ e curvas de poder

### Notas de calibração
- O null da lei de Benford **não descreve contagens eleitorais limpas**
  (erro tipo I até 0,72 com seções homogêneas). O default `null = :resampled`
  calibra o p-valor contra a distribuição empírica de contagens;
  `null = :benford` fica disponível como estatística descritiva.
- O teste do **penúltimo dígito** é opt-in: seu null uniforme rejeita eleições
  limpas em até 100 % dos casos quando as contagens cobrem menos de duas
  décadas. `null_valid` sinaliza a condição.
- Nenhum teste da bateria detecta *ballot stuffing* proporcional.
```

---

## 9. Como reexecutar

```bash
julia --project=validation -e 'using Pkg
  Pkg.develop(path="."); Pkg.add(["Distributions","HypothesisTests","Aqua","JET"])'

julia --project=validation -t auto validation/01_calibracao_benford.jl
julia --project=validation -t auto validation/07_penultimo_digito.jl   # C1
julia --project=validation -t auto validation/13_benford_null_reamostrado.jl  # C2
python3 -m venv .venv && .venv/bin/pip install scipy numpy
julia --project=validation validation/10a_export_para_scipy.jl
.venv/bin/python validation/10b_xref_scipy.py                          # referência cruzada
```

O índice completo dos scripts está em `validation/README.md`.
