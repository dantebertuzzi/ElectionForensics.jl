# Revisão técnica — ElectionForensics.jl

**Revisor:** auditoria de engenharia de pacotes Julia + estatística forense eleitoral
**Commit auditado:** `590903f` · **Correções:** `051d375` (branch `revisao-correcoes`)
**Data:** 2026-08-29 · **Julia:** 1.12.7 (também verificado em 1.10 LTS e 1.6.7)
**Scripts de evidência:** `validation/` (16 scripts reexecutáveis)

---

## 1. Veredito

> ## Antes: `NÃO APTO` · Depois das correções: `APTO COM RESSALVAS`

**Todas as correções foram aplicadas e commitadas** na branch
`revisao-correcoes`. Dos 26 achados, 25 estão fechados e 1 continua parcial
(M3). O que segue descreve o pacote **como auditado** (`590903f`); o estado
atual está na seção 5.

### O diagnóstico original

Dois dos quatro testes da bateria são **anti-conservadores por construção**: no
DGP eleitoral mais comum no Brasil (seções de tamanho homogêneo), o teste do
penúltimo dígito rejeita eleições limpas em **50–100 %** dos casos e o teste de
Benford 2BL em **72 %** — e ambos são executados por default em
`forensics_report`. Um pacote de forense eleitoral que acusa fraude em dados
limpos com essa frequência não pode ser publicado no General Registry.

A boa notícia: **a engenharia estava correta e as correções eram localizadas.**
A maquinaria χ², as probabilidades de Benford, os limiares de Nigrini, a fórmula
do p-valor Monte Carlo e o teste de Rozenas conferem com a literatura e com
implementações independentes.

### As ressalvas que permanecem

1. **M3 está fechado** — mas não por um default melhor. Nenhum `min_value`
   calibra todos os regimes e nenhum diagnóstico escalar prediz a falha
   (`validation/18`: correlações de +0,25 a −0,55, com contraexemplos nos dois
   sentidos). A resposta é `calibration_check`: o pacote **mede** a taxa de
   erro tipo I nos dados do usuário e avisa quando o p-valor não é confiável.
   Ver seção 5c.
2. **A validação em dados reais usa fontes estrangeiras.** O CDN do TSE
   bloqueia clientes não-browser deste ambiente, então a calibração real foi
   feita em 29 pleitos americanos (OpenElections) e o poder na Rússia 2012
   (Shpilkin / Kobak et al.). Ver seção 5b. O TSE segue não verificado.
3. **A limitação de exatidão do teste de Rozenas permanece** — é inerente ao
   método (ver 5b), está documentada, e `max_denom` maior a mitiga.
4. **As correções são minhas, não revisadas por terceiro.** Foram verificadas
   por execução (178 testes, Aqua, JET, Documenter, Julia 1.10 e 1.12), mas
   passaram por um único par de olhos.

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

## 5. Estado após as correções

Commit `051d375` na branch `revisao-correcoes`. Verificado por execução:

```
Test Summary:        | Pass  Total   Time
ElectionForensics.jl |  178    178  35.1s     (era 76)
     Testing ElectionForensics tests passed
```

`Aqua.test_all` sem exclusões · JET sem relatórios em código do pacote ·
Documenter compila com `checkdocs = :exports` e `doctest = true` ·
178/178 também em Julia 1.10 LTS, do zero.

### Calibração do código entregue (α = 0,05, m = 2000)

| Cenário (eleições **limpas**) | 2BL antes | 2BL depois | Últ. dígito | Rozenas |
|---|---|---|---|---|
| totals 200–400 (caso TSE) | **0,710** | **0,054** | 0,050 | 0,050 |
| totals 150–900 | 0,332 | 0,038 | 0,052 | 0,048 |
| shares Beta(0.5, 0.5) | — | 0,064 | 0,052 | 0,052 |
| shares Beta(0.3, 0.3) | — | 0,050 | **0,082** ⚠️ | 0,030 |

O penúltimo dígito saiu do default: era 0,503 e 1,000 nos dois primeiros
cenários.

### Poder do teste de Rozenas (metas redondas)

| ε | Beta(8,6) antes | Beta(8,6) depois | Beta(0.5,0.5) antes | Beta(0.5,0.5) depois |
|---|---|---|---|---|
| 0,00 | 0,040 | 0,034 | 0,036 | 0,036 |
| 0,01 | 0,196 | 0,188 | 0,168 | **0,288** |
| 0,02 | 0,472 | **0,540** | 0,480 | **0,700** |
| 0,05 | 0,988 | 0,988 | 0,986 | **1,000** |

### M5: minha hipótese original estava errada

No relatório inicial atribuí o conservadorismo do teste de Rozenas à **dupla
adição de ruído binomial** e propus a correção de variância de Silverman &
Young (1987). O experimento (`validation/15_rozenas_escala.jl`) refutou isso:

| dist. de shares | atual | `shrink` (minha proposta) | banda menor | **logit** |
|---|---|---|---|---|
| Beta(0.5, 0.5) | 0,012 | 0,012 | 0,032 | **0,036** |
| Beta(0.3, 0.3) | 0,000 | 0,002 | 0,016 | **0,030** |
| poder ε = 2 % | 0,480 | 0,424 | 0,638 | **0,666** |

A correção de variância não move nada — o desvio que ela corrige é de ~3 %. A
causa real é a **banda global vazando nas fronteiras**: com `h ≈ 0,07` sob
Beta(0.3,0.3), o jitter empurra seções quase unânimes para o miolo, onde as
frações coarse são densas. Na escala logit o passo encolhe sozinho perto de 0 e
1. A correção implementada é essa, não a que eu havia proposto.

### Achados fechados nesta rodada

| Sev. | # | O que foi feito |
|---|---|---|
| 🔴 | C1 | `:penultimate` opt-in, campo `null_valid`, aviso e docstring corrigida |
| 🔴 | C2 | `null = :resampled` como default; `:benford` mantido como descritivo |
| 🟠 | A1 | `Manifest.toml` desversionado + `.gitignore` |
| 🟠 | A2 | `[compat]` para `Printf`, `Random`, `Statistics`, `Tables`, `Test` |
| 🟠 | A3 | Exemplo do README corrigido |
| 🟠 | A4 | Campo `qvalues` (Benjamini–Hochberg) exibido ao lado do z |
| 🟠 | A5 | `clamp` → reflexão, depois superado pela escala logit |
| 🟠 | A6 | CI (3 SOs × 3 versões), TagBot, CompatHelper, `docs/` que compila |
| 🟠 | A7 | Docstrings nos três tipos exportados |
| 🟡 | M1 | Regra de Cochran (`min_expected`) com aviso |
| 🟡 | M2 | Uma semente por testset; asserções sobre taxa; valores analíticos |
| 🟡 | M4 | Coerção explícita: `Float64` inteiro ok, fracionário e `NaN` erram, `missing` exige `skipmissing` |
| 🟡 | M5 | Jitter na escala logit (`boundary = :logit`) |
| 🟡 | M6 | Aviso para `m < 30` |
| 🟡 | M7 | `_H_FLOOR` nomeado; `_silverman` devolve `(h, degenerate)` |
| 🟡 | M8 | README qualifica a aproximação e remete aos `qvalues` |
| 🟡 | M9 | Interface Tables.jl nos três resultados, com `schema` explícito |
| 🟡 | M10 | `forensics_report` pula testes inaplicáveis em vez de abortar |
| 🟡 | M11 | O `show` explica por que MAD e χ² podem discordar |
| 🟢 | B1 | Fronteiras de Nigrini fechadas no limite superior |
| 🟢 | B2 | `_digit_counts!` parametrizado — sem `Union` de tipos de função |
| 🟢 | B3 | PNGs desversionados |
| 🟢 | B4 | `warn` como keyword |
| 🟢 | B5 | Chave `(num, den)` no lugar de `hash(::Rational)` no laço quente |
| 🟢 | B6 | `min_value` abaixo do mínimo estrutural agora é erro |

### O que continua aberto

| Sev. | # | Situação |
|---|---|---|
| 🟡 | M3 | **Parcial.** Aviso de dispersão e guarda de `min_value` adicionados, mas o teste do último dígito ainda roda e chega a 0,082 sob shares bimodais. Fechar exigiria um null reamostrado próprio ou um filtro por `sd`. |
| — | — | Dados reais do TSE: `HTTP 403` deste ambiente. Script pronto em `validation/14_dados_reais_tse.jl`. |

---

## 5b. Validação em dados reais

Duas fontes públicas de proveniência documentada, ambas reproduzíveis por
`validation/fetch_dados_reais.py`:

| Fonte | O que é | Testa |
|---|---|---|
| **OpenElections** | Resultados por precinct das gerais americanas de 2020, transcritos dos boletins oficiais das secretarias eleitorais estaduais | Calibração |
| **dkobak/elections** | 95 413 seções (UIK) da presidencial russa de 2012, coletadas por Sergey Shpilkin — é o dataset de Kobak, Shpilkin & Pshenichnikov (2016), *Annals of Applied Statistics* 10(1):54–73 | Poder, contra verdade publicada |

### Calibração: 29 pleitos americanos

29 dos 35 estados baixados passaram no crivo de qualidade da extração; os
outros 6 foram excluídos e o motivo está registrado em `qualidade.csv`.

| Teste | rej @ α=0,05 | rej @ α=0,10 | KS vs U(0,1) |
|---|---|---|---|
| 2BL, **null clássico** (lei de Benford) | **8/29 — 28 %** | 12/29 — 41 % | **0,0001** |
| 2BL, **null reamostrado** (novo default) | **0/29 — 0 %** | 2/29 — 7 % | 0,160 |
| Último dígito | 3/29 — 10 % | 5/29 — 17 % | 0,273 |
| Frações coarse (Rozenas) | 3/29 — 10 % | 3/29 — 10 % | 0,420 |

**Este é o teste mais forte do relatório, e confirma o achado C2 em dados
reais.** O null da lei de Benford rejeita 28 % de eleições americanas
legítimas — NY com p = 3,7 × 10⁻¹⁹, MA com 7,8 × 10⁻¹¹, IL com 2,0 × 10⁻⁹,
OH com 1,5 × 10⁻⁸ — e a distribuição dos p-valores é incompatível com a
uniformidade. O null reamostrado não rejeita nenhum, com KS não significativo.
A simulação previa 0,71 → 0,054; os dados reais entregam 0,28 → 0,00.

Duas ressalvas honestas sobre a mesma tabela:

- O teste do último dígito dá 3/29 contra 1,5 esperado (ID, LA, MO). Não é
  significativo (KS = 0,273), mas é consistente com o resíduo do achado M3.
- As frações coarse dão 3/29 (AK, SD, LA). Eu havia levantado, com uma amostra
  de 8 pleitos, a hipótese de **sobredispersão da nula em dados reais**.
  Com 29 pleitos ela **não se sustenta**: KS = 0,420. Era ruído amostral.

### Poder: Rússia 2012

| Teste | Resultado |
|---|---|
| Benford 2BL clássico | p = 4,8 × 10⁻¹¹ |
| Benford 2BL reamostrado | p = 0,012 |
| Último dígito | p = 0,059 |
| **Frações coarse** | **p = 0,043** · T = 2122 · E[T\|H₀] = 2051 · z = 1,7 |

As duas frações com maior excesso são **4/5 (z = 2,9)** e **3/4 (z = 2,6)** —
80 % e 75 %, exatamente os percentuais redondos que Kobak et al. documentam.
O pacote detecta a anomalia canônica da literatura, mas fracamente.

**Por que fracamente — uma limitação inerente que os dados reais expuseram.**
O primitivo exige que o percentual seja *exatamente* k/d, o que requer
`totals` divisível por d. Nos dados russos:

| d | seções com d \| N | % |
|---|---|---|
| 2 | 47 807 | 50,1 |
| 3 | 31 843 | 33,4 |
| 5 | 19 705 | 20,7 |
| 10 | 10 053 | 10,5 |

77,2 % das seções são elegíveis para *alguma* fração com d ≤ 10, mas só 2,6 %
acertam uma. Kobak et al. comparam a densidade em percentuais inteiros contra a
vizinhança local e não precisam de divisibilidade — por isso enxergam um efeito
muito maior nos mesmos dados. Aumentar `max_denom` recupera parte do sinal
(p = 0,042 com d ≤ 10 → 0,017 com d ≤ 100), o que sugere revisar o default.

### Confirmação de C1 em dados reais

O teste do penúltimo dígito devolve **p = 3,4 × 10⁻⁶³** na Rússia 2012 — e
`null_valid = false`, porque as contagens cobrem 1,31 décadas. O **mesmo** teste
devolve p ≈ 10⁻⁸ em dados simulados *limpos*. O número não distingue fraude de
artefato do null; agora é opt-in e vem marcado.

### Referência cruzada em dados reais

16 comparações de χ², p-valor e MAD contra `scipy` sobre os mesmos datasets
reais: discrepância relativa máxima **2,8 × 10⁻¹²**. MAD contra `benford_py`
(biblioteca Python independente): **1 × 10⁻¹⁴**.

### Três erros meus na preparação dos dados

Registrados porque são o tipo de coisa que produz "fraude" inexistente:

1. Uma regex removeu **Don Blankenship**, candidato real de 2020, por o
   sobrenome conter "blank".
2. Linhas de subtotal por condado (`"ANDROSCOGGIN Total"`) entraram como
   precinct, duplicando seções.
3. A coluna `absentee_votes` de Rhode Island foi ignorada, encolhendo o
   denominador.

Os três alteravam resultados. Por isso `fetch_dados_reais.py` agora verifica,
por estado, a identidade `linha-resumo = Σcandidatos + brancos/nulos` e a
coerência do top-2, e exclui quem não passa.


---

## 5c. Como a ressalva M3 foi resolvida

### O defeito, demonstrado em dados reais

Dos 29 pleitos americanos, três rejeitaram no teste do último dígito. A
Louisiana é a demonstração limpa do problema:

| Corte | p-valor | Perfil dos dígitos 0→9 (%) |
|---|---|---|
| `min_value = 10` (default) | **0,0034** | 11,6 · 10,1 · 10,8 · 10,5 · 10,2 · 10,2 · 9,3 · 9,6 · 8,6 · 9,1 |
| `min_value = 100` | **0,479** | 11,1 · 9,8 · 10,4 · 10,7 · 9,8 · 10,1 · 9,8 · 9,5 · 8,9 · 10,0 |

O perfil é um **gradiente monótono decrescente**, não um excesso em 0 e 5 — a
assinatura de artefato do null, não de fabricação humana. 20,8 % das contagens
retidas da Louisiana estão abaixo de 50. Restringir a contagens ≥ 100 dissolve
a rejeição.

### Por que não bastava mudar o default

| Regime (eleições limpas) | mv=10 | mv=20 | mv=30 | mv=50 | mv=100 | null MC |
|---|---|---|---|---|---|---|
| Beta(8,6) totais 150:900 | 0,050 | 0,050 | 0,050 | 0,049 | 0,053 | 0,048 |
| Beta(0.3,0.3) totais 150:900 | 0,078 | 0,059 | 0,065 | 0,059 | 0,044 | 0,072 |
| Beta(8,6) totais 20:400 | 0,052 | 0,075 | 0,073 | 0,073 | **0,077** | 0,048 |
| Beta(0.3,0.3) totais 20:400 | **0,136** | 0,099 | 0,081 | 0,062 | 0,049 | 0,111 |

`min_value = 100` conserta três regimes e **piora** um. Um null reamostrado —
que funcionou para a 2BL — não ajuda aqui (0,111 no pior caso). E não há
diagnóstico escalar que prediga a falha: testei proporção de contagens
pequenas, 5º percentil e desvio-padrão relativo, com correlações de Spearman de
+0,245, −0,545 e −0,027 contra o erro tipo I, e contraexemplos em ambas as
direções (`validation/18_criterio_ultimo_digito.jl`).

### A resposta: o pacote mede em vez de prometer

`calibration_check(votes, totals)` simula eleições limpas com os **seus**
`totals` e um formato de percentuais como o seu, roda cada teste, e reporta a
taxa de rejeição que ele de fato entrega.

```
Auto-diagnóstico de calibração
2500 seções · 200 eleições limpas simuladas · α = 0.050
  teste            rejeição   esperado    veredito
  Benford          0.027      0.050±0.031 calibrado
  último dígito    0.367      0.050±0.031 anticonservador
  penúltimo dígito 0.960      0.050±0.031 anticonservador
  frações coarse   0.080      0.050±0.031 calibrado
⚠ testes anticonservadores rejeitam dados LIMPOS acima de α
```

**O diagnóstico acerta** (`validation/19_valida_calibration_check.jl`):

| Regime | Erro tipo I medido | Previsto pelo diagnóstico | Veredito |
|---|---|---|---|
| Beta(8,6) totais 150:900 | 0,040 | 0,060 | calibrado |
| Beta(0.3,0.3) totais 150:900 | 0,068 | 0,060 | calibrado |
| Beta(0.3,0.3) totais 20:400 | 0,100 | 0,136 | anticonservador |
| Beta(8,6) totais 10:120 | **0,372** | **0,364** | anticonservador |
| Beta(8,6) totais 1000:9000 | 0,052 | 0,028 | calibrado |

**E separa corretamente os três casos reais:**

| UF | p (último dígito) | Taxa medida em dados limpos como os dela | Veredito |
|---|---|---|---|
| **LA** | 0,0034 | **0,164** | anticonservador — p-valor não confiável |
| ID | 0,0203 | 0,068 | calibrado — p-valor confiável |
| MO | 0,0215 | 0,064 | calibrado — p-valor confiável |
| NY | 0,3194 | 0,068 | calibrado |
| WY | 0,7189 | 0,032 | calibrado |

A Louisiana — o falso positivo que eu havia demonstrado à mão — é marcada.
Idaho e Missouri passam, e são exatamente o ruído de teste múltiplo esperado
(1,5 rejeições esperadas em 29 a α = 0,05; sobram 2).

`forensics_report` roda a checagem por default e a imprime **antes** de
qualquer p-valor.


---

## 5d. Comparação com o estado da arte

### Em Julia não há concorrente

Varredura do General Registry: **14 249 pacotes, nenhum de forense eleitoral ou
de lei de Benford**. Os únicos nomes que casam com termos do domínio são
`SynchronicBallot` (consenso distribuído) e `BRElections` (acesso a dados do
TSE, do mesmo autor). `ElectionForensics.jl` seria o primeiro.

### As bibliotecas de Benford têm a mesma ressalva — e pior

Rodei cada biblioteca **com seus próprios defaults** sobre os mesmos 29 pleitos
americanos legítimos (`validation/20_concorrentes.py`):

| Biblioteca | Acusa | Taxa | KS vs U(0,1) |
|---|---|---|---|
| `benford.analysis` (CRAN) — MAD 2 díg. | **28/29** | **97 %** | sem p-valor |
| `BenfordTests::chisq.benftest` (1BL) | 23/29 | 79 % | 3 × 10⁻¹⁸ |
| `benfordslaw` (PyPI, default) | 23/29 | 79 % | 3 × 10⁻¹⁸ |
| `benford_py` (PyPI) — MAD 1BL | 23/29 | 79 % | sem p-valor |
| `BenfordTests::meandigit.benftest` | 21/29 | 72 % | 4 × 10⁻¹⁴ |
| `BenfordTests::ks.benftest` | 20/29 | 69 % | 1 × 10⁻¹⁶ |
| `BenfordTests::chisq.benftest` (2 díg.) | 19/29 | 66 % | 1 × 10⁻¹² |
| **ElectionForensics.jl** — null clássico | 8/29 | 28 % | 1 × 10⁻⁴ |
| **ElectionForensics.jl** — null reamostrado | **0/29** | **0 %** | 0,16 |

`benford_py` e `benfordslaw` foram **executados**. As estatísticas de
`BenfordTests` e `benford.analysis` foram **reimplementadas a partir do fonte R**
(não há R neste ambiente) — os limiares de MAD conferem com
`internal.functions-new-code-2.R` (`0.0012, 0.0018, 0.0022`).

Duas observações de justiça:

- Essas são ferramentas **de propósito geral** para auditoria contábil, não de
  forense eleitoral. Aplicá-las a contagens de votos é uso indevido — mas é um
  uso que elas convidam: nenhuma verifica a condição de aplicabilidade.
  `benfordslaw` é a única que ao menos **documenta** a premissa ("the numbers
  should cover several orders of magnitude", docstring), sem checá-la.
- O pacote já saía na frente por um acerto do autor: o default `digit = 2`
  (2BL, seguindo Mebane 2008) rende 28 % contra os 79 % da 1BL. A escolha
  original estava certa; só não era suficiente.

### `spikes` (Rozenas) é melhor onde importa — e tem outras ressalvas

`spikes` v1.1 (CRAN, 2016) é a implementação do **próprio autor do método**.
Lendo o fonte:

| | `spikes` | `ElectionForensics.jl` |
|---|---|---|
| Densidade latente | Mistura de Beta-Binomiais ajustada por **EM** | Kernel gaussiano sobre os shares observados |
| Reamostragem | Posterior `Beta(y+α, n−y+β)` — **deconvolui o ruído binomial** | Jitter em escala logit |
| Fronteiras 0/1 | Naturais na Beta | Precisou de tratamento explícito (achado M5) |
| Saída | % estimada de seções fraudulentas + IC bootstrap | p-valor Monte Carlo + q-valores BH |
| Controle de erro tipo I | **Nenhum** (não é teste) | Medido por `calibration_check` |
| Escopo | Só frações coarse | + Benford, último/penúltimo dígito |
| Manutenção | Última versão 2016 | — |
| Efeitos colaterais | `plot()` dentro do ajuste, `dev.off()` incondicional, laços `while` sem teto | — |

**O modelo de nula de `spikes` é superior ao deste pacote** — é exatamente a
deconvolução que meu experimento do M5 tateava. Recomendo, como trabalho
futuro, substituir o jitter por uma mistura Beta-Binomial ajustada aos dados.

`eforensics` (implementação do modelo de Mebane & Klašnja) existe no GitHub mas
nunca foi lançado: o `DESCRIPTION` ainda traz o texto-modelo do roxygen
("What the Package Does (one line, title case)"), versão 0.0.0.9000, fora do CRAN.

### O diferencial

Nenhuma das bibliotecas examinadas — nem as de Benford, nem `spikes` — mede a
própria taxa de erro tipo I nos dados do usuário. `calibration_check` não é um
teste melhor; é o único que avisa quando o próprio p-valor não é confiável.


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
- [x] `[compat]` para `Printf`, `Random`, `Statistics`, `Test` *(patch)*
- [x] `[extras]` / `[targets]` de teste declarados
- [n/a] `[weakdeps]` / `[extensions]` — o pacote não tem extensões
- [x] `julia = "1.6"` **verificado** (76/76 em 1.6.7) — o patch sobe para `"1.10"` (LTS atual)
- [x] `Manifest.toml` fora do controle de versão *(patch)*
- [x] `.gitignore` *(patch)*

**Documentação**
- [x] `LICENSE` presente (MIT)
- [x] README com seção de métodos e referências
- [x] README com exemplo **executável** *(patch: `max_denom` corrigido)*
- [x] Docstrings nas 9 funções exportadas
- [x] Docstrings nos 3 tipos exportados *(patch)*
- [x] `docs/` com Documenter — criado, compila com `checkdocs = :exports`
- [x] Doctests — `docs/src/index.md` tem `jldoctest`, rodam no build e na CI
- [x] CHANGELOG

**CI**
- [x] Matriz Julia (1.10 LTS, 1, pre) *(patch)*
- [x] Múltiplos SOs (Linux, macOS, Windows) *(patch)*
- [x] TagBot, CompatHelper *(patch)*
- [x] Cobertura (Codecov) *(patch)*
- [x] Job de doctests com `docs/` presente

**Qualidade**
- [x] `Aqua`: ambiguidades, pirataria, deps obsoletas, exports — limpo
- [x] `Aqua`: `[compat]` *(patch)*
- [x] `JET`: nada em código do pacote
- [x] `] test` passa do zero em ambiente limpo (1.6.7, 1.10.12, 1.12.7)
- [x] Precompila sem warnings
- [x] Tempo de carga aceitável (0,38 s)
- [x] `@warn` suprimível *(patch: keyword `warn`)*

**Correção estatística**
- [x] Probabilidades de Benford conferidas contra a literatura
- [x] Graus de liberdade corretos
- [x] Referência cruzada com `scipy` (Δ ≤ 2,8 × 10⁻¹²)
- [x] p-valor Monte Carlo com a convenção correta
- [x] RNG reprodutível, sem dependência de threads
- [x] Calibração sob H₀ da 2BL *(patch: 0,710 → 0,058)*
- [x] Penúltimo dígito fora do default *(patch)*
- [x] Multiplicidade nas frações *(patch: BH)*
- [x] Conservadorismo do Rozenas sob shares em U — corrigido pela escala logit (0,000 → 0,030; poder 0,48 → 0,70)
- [x] Retornos compatíveis com Tables.jl — os três tipos, com `schema` explícito

- [x] **M3 resolvido** por `calibration_check` (seção 5c) — o pacote mede a taxa de erro tipo I nos dados do usuário e marca os testes não confiáveis

**Não verificado**
- [ ] Dados reais do TSE — o CDN respondeu **HTTP 403** a clientes não-browser
      neste ambiente (`cdn.tse.jus.br` e `dadosabertos.tse.jus.br`). Script
      pronto em `validation/14_dados_reais_tse.jl`; rode localmente.
- [x] Build do Documenter — compila localmente sem erro nem warning de docstring.

**Não aplicável (não existe no pacote)**
- Modelo de Klimek et al. / fingerprint 2D · quadratura de Gauss–Hermite ·
  resampling bayesiano de Rozenas · `logsumexp` · bootstrap percentil/BCa ·
  paralelização por threads · extensões Makie · tema Dracula

---

## 8. Versão sugerida e changelog

O pacote **nunca foi publicado** (uma única tag e nenhum release no General).
Não há bump a fazer: `0.1.0` continua correto, agora já com as correções.
Lançar 0.1.0 com os defeitos e corrigir em 0.2.0 teria colocado no registro
permanente uma versão que acusa fraude em dados limpos.

`CHANGELOG.md` foi criado no repositório com o conteúdo abaixo.

Veja `CHANGELOG.md`.

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
