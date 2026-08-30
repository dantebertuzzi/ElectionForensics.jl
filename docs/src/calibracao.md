# Calibração

```@meta
CurrentModule = ElectionForensics
```

Um teste de fraude que rejeita eleições limpas é pior que nenhum teste. Esta
página resume o que foi medido por simulação; os scripts estão em `validation/`.

## A lei de Benford não descreve contagens eleitorais

Sob eleições limpas simuladas (`votos ~ Binomial(N, p)`, `p ~ Beta`), o χ²
contra a lei de Benford rejeita H₀ com frequência muito acima de α:

| Tamanho das seções | 1BL | 2BL clássica | 2BL reamostrada |
|---|---|---|---|
| 200–400 (típico do TSE) | 1,000 | **0,710** | 0,058 |
| 150–900 | 1,000 | 0,332 | 0,038 |
| log-normal(6; 1,2) | 0,057 | 0,062 | 0,048 |

A conformidade com Benford só aparece quando os tamanhos de seção são
aproximadamente log-normais com dispersão grande — a crítica de Deckert,
Myagkov & Ordeshook (2011). A maquinaria χ² do pacote está correta: com dígitos
sorteados da própria lei de Benford, a taxa de rejeição fica em 0,046–0,054.
O problema é o modelo, não o código.

Por isso o default é [`BenfordConfig`](@ref) com `null = :resampled`, que
compara o dígito observado contra a **distribuição empírica de contagens
suavizada em `log10`**. `null = :benford` continua disponível como estatística
descritiva (MAD e conformidade de Nigrini).

## O penúltimo dígito não é uniforme

A uniformidade do último dígito exige densidade suave numa escala de 10 — vale
na prática. O penúltimo exige suavidade numa escala de 100, o que quase nunca
vale:

| Cenário (eleições limpas) | Erro tipo I do penúltimo |
|---|---|
| totais 150–900 | 0,503 |
| totais 100–200 | **1,000** |
| totais 100–9999 | 0,053 |

Sob H₀ a distribuição do penúltimo dígito é **decrescente** (10,95 % → 8,87 %),
não uniforme. Um null reamostrado reduz mas não elimina a inflação. Por isso
`:penultimate` é opt-in em [`forensics_report`](@ref) e o resultado carrega
`null_valid`.

O teste do **último** dígito, em contraste, é bem calibrado (0,044–0,058 para
seções de 60 a 5 000 eleitores) e tem poder 1,000 contra arredondamento decimal
em 5 % das seções.

## A escala do jitter no teste de Rozenas

A nula reamostra percentuais com um kernel gaussiano. Com a banda global de
Silverman na escala de probabilidade, seções quase unânimes têm sua massa
empurrada para o miolo, onde as frações coarse são densas — a nula superestima
`T` e o teste fica conservador. O jitter na escala **logit** encolhe sozinho
perto das fronteiras:

| Shares | erro tipo I `:reflect` | erro tipo I `:logit` | poder (ε = 2 %) |
|---|---|---|---|
| Beta(8,6) | 0,032 | 0,038 | 0,456 → 0,474 |
| Beta(0.5, 0.5) | 0,012 | 0,036 | 0,480 → **0,666** |
| Beta(0.3, 0.3) | 0,000 | 0,030 | — |

`boundary = :logit` é o default; `:reflect` preserva o comportamento anterior.

## Meça, não confie: `calibration_check`

As três seções acima mostram que a validade de cada null depende do **formato**
da distribuição de contagens, e que nenhuma estatística escalar simples prediz
a falha — `validation/18_criterio_ultimo_digito.jl` testa três candidatas e
encontra correlações de +0,25 a −0,55, com contraexemplos nos dois sentidos.

Por isso o pacote não promete calibração: ele a **mede** nos seus dados.
[`calibration_check`](@ref) simula eleições limpas com os seus tamanhos de
seção e um formato de percentuais como o seu, roda cada teste, e reporta a taxa
de rejeição que ele de fato entrega.

```julia
julia> calibration_check(votes, totals)
Auto-diagnóstico de calibração
2500 seções · 200 eleições limpas simuladas · α = 0.050
  teste            rejeição   esperado    veredito
  Benford          0.027      0.050±0.031 calibrado
  último dígito    0.367      0.050±0.031 anticonservador
  penúltimo dígito 0.960      0.050±0.031 anticonservador
  frações coarse   0.080      0.050±0.031 calibrado
⚠ testes anticonservadores rejeitam dados LIMPOS acima de α; não conclua
  fraude a partir deles nestes dados
```

O diagnóstico acerta: em `validation/19_valida_calibration_check.jl`, o pior
regime tem erro tipo I **medido** de 0,372 e **previsto** de 0,364.

[`forensics_report`](@ref) roda essa checagem por default (`calibrate = true`)
para os testes de dígito e imprime o resultado **antes** de qualquer p-valor.

## Duas nulas para o teste de Rozenas

`null = :kernel` (default) perturba os percentuais observados com um kernel
gaussiano em escala logit. `null = :betabinomial` implementa o modelo do pacote
`spikes`, do próprio Rozenas: ajusta por EM uma mistura de Beta-Binomiais à
distribuição **latente** de percentuais e reamostra da posterior de cada seção,

```math
\tilde p_i \sim \mathrm{Beta}(y_i + \alpha_k,\; n_i - y_i + \beta_k), \qquad
y_i^* \sim \mathrm{Binomial}(n_i, \tilde p_i)
```

com o componente `k` sorteado proporcionalmente à responsabilidade. A posterior
deconvolui o ruído binomial e vive em [0,1], então dispensa tratamento de
fronteira.

**Em simulação as duas são equivalentes** — mesmo quando o processo gerador
é, por construção, uma mistura de Betas:

| Shares | erro tipo I `:kernel` | `:betabinomial` | poder ε=2 % `:kernel` | `:betabinomial` |
|---|---|---|---|---|
| Beta(8,6) | 0,030 | 0,033 | 0,420 | 0,443 |
| Beta(0.5,0.5) | 0,038 | 0,033 | 0,515 | 0,530 |
| bimodal 2 comp. | 0,033 | 0,045 | — | — |

**Em dados reais a versão paramétrica pode falhar.** Nos 95 413 UIK russos de
2012 a mistura ajusta o corpo da distribuição muito bem (razão observado/predito
entre 0,97 e 1,05 para percentuais acima de 0,3) mas subestima a cauda esquerda
em **3,8×** — e a nula passa a "encontrar" excesso justamente nas frações que
caem ali: 1/7, 1/6, 1/5 e 1/4 aparecem com q < 0,05 sem que haja nada lá.
A nula por kernel, sendo não-paramétrica, acompanha a distribuição observada
inclusive nas caudas e não tem esse modo de falha.

Por isso `:kernel` continua o default, e toda mistura ajustada carrega o
diagnóstico `misfit`: a maior razão observado/predito entre as faixas de
percentual com massa não desprezível. Acima de 2 o pacote avisa. Nos dados
russos ele marca 2,57 na faixa [0,20, 0,25) — exatamente onde estão 1/5 e 1/4.

## O fingerprint não é um teste

[`election_fingerprint`](@ref) não produz p-valor, e assimetria/curtose
isoladas não são teste de hipótese: com dezenas de milhares de seções qualquer
desvio ínfimo da normalidade rende uma estatística enorme. Klimek et al. as
usam como descrição comparativa entre países, e é assim que devem ser lidas.

Reprodução da Fig. 3 do artigo com dados reais
(`validation/22_fingerprint.jl`):

| Pleito | Assimetria de ν | Curtose | Canto ≥95 %/≥95 % |
|---|---|---|---|
| simulada limpa (20 000 seções) | +0,13 | 3,22 | 0 |
| EUA 2020 · Illinois | **−0,01** | **3,00** | 0 |
| EUA 2020 · Alabama | +1,53 | 4,78 | 0 |
| Rússia 2012 | **−2,15** | **9,12** | **2 677 (2,81 %)** |

Illinois cai em cima do (0, 3) previsto. Alabama se afasta sem nenhuma seção no
canto — o artigo atribui desvios assim a heterogeneidade da população (o
exemplo dele é o Canadá), e aqui há ainda viés de subamostra, porque só 1 204
dos 2 110 precincts publicam eleitorado. A Rússia mostra as duas assinaturas
que o artigo descreve.

## Multiplicidade

`RozenasResult` reporta excessos por fração. Sob H₀, o maior z entre as 31
frações tem mediana 2,41, e `P(algum z > 2) = 0,678`. Use o campo `qvalues`
(Benjamini–Hochberg), nunca os z brutos.
