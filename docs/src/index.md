# ElectionForensics.jl

Forense eleitoral em Julia: lei de **Benford** (1º e 2º dígitos), teste de
**último dígito** (Beber & Scacco 2012) e teste de **frações coarse**
(Rozenas 2017), com relatório integrado no terminal.

```@meta
CurrentModule = ElectionForensics
```

```@docs
ElectionForensics
```

## Instalação

```julia
] add https://github.com/dantebertuzzi/ElectionForensics.jl
```

## Uso

```jldoctest exemplo
julia> using ElectionForensics, Random, Distributions

julia> rng = MersenneTwister(2026);

julia> totals = rand(rng, 200:400, 1500);

julia> p = rand(rng, Beta(8, 6), 1500);

julia> votes = [rand(rng, Binomial(totals[i], p[i])) for i in 1:1500];

julia> r = benford_test(votes; digit = 2, B = 199, rng = rng, warn = false);

julia> r.digit, r.n, length(r.counts)
(2, 1500, 10)

julia> sum(r.counts) == r.n
true
```

Os resultados são tabelas: qualquer consumidor de Tables.jl os aceita.

```jldoctest exemplo
julia> using Tables

julia> Tables.schema(r).names
(:digit, :count, :observed, :expected, :excess, :pearson_residual)

julia> length(Tables.rows(r))
10
```

Para a bateria completa:

```julia
res = forensics_report(votes, totals)
res.rozenas.pvalue
res.benford.conformity
```

## O que estes testes não provam

Nenhum teste isolado prova fraude, e **nenhum deles detecta *ballot stuffing*
proporcional**: em simulação, adicionar 15 % dos votos a 20 % das seções deixa
os quatro testes exatamente no nível nominal. O teste de Rozenas detecta metas
percentuais redondas; o de último dígito detecta arredondamento decimal. Fraude
sem assinatura em dígitos ou em frações exatas passa despercebida.

Leia também a página [Calibração](calibracao.md): dois dos nulls clássicos da
literatura **não** são válidos para contagens eleitorais, e o pacote se afasta
deles deliberadamente.
