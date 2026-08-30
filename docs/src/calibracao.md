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

## Multiplicidade

`RozenasResult` reporta excessos por fração. Sob H₀, o maior z entre as 31
frações tem mediana 2,41, e `P(algum z > 2) = 0,678`. Use o campo `qvalues`
(Benjamini–Hochberg), nunca os z brutos.
