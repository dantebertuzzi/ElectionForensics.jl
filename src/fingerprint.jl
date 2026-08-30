# ── fingerprint eleitoral (Klimek et al. 2012) ───────────────────────
#
# Klimek, P., Yegorov, Y., Hanel, R. & Thurner, S. (2012). "Statistical
# detection of systematic election irregularities". PNAS 109(41):16469–16473
# (preprint arXiv:1201.3087).
#
# Duas coisas do artigo são implementadas aqui, ambas especificadas por
# completo no texto principal:
#
#  1. O **fingerprint**: histograma 2-D do número de seções por comparecimento
#     (eixo x) e por percentual de votos no vencedor (eixo y). Em eleições
#     limpas há um único aglomerado; manipulação incremental o "borra" em
#     direção ao canto superior direito, e manipulação extrema cria um segundo
#     aglomerado junto de (100 %, 100 %).
#
#  2. A **taxa logarítmica de voto** νᵢ = log((Nᵢ − Wᵢ)/Wᵢ), cuja distribuição
#     reescalada é aproximadamente gaussiana em eleições limpas. Assimetria e
#     curtose de ν são o diagnóstico da Fig. 3 do artigo: os países sem
#     alegação de fraude se agrupam em torno de (0, 3), Rússia e Uganda se
#     afastam muito.
#
# O modelo paramétrico de fraude (fᵢ, f_e) do artigo NÃO está implementado: a
# especificação completa dos mecanismos está no Supporting Information, que não
# acompanha o preprint. Implementá-lo de memória seria inventar o modelo.

"""
    FingerprintConfig

Parâmetros do [`election_fingerprint`](@ref).

- `bins`: resolução do histograma 2-D em cada eixo.
- `vote_axis`: denominador do eixo vertical.
  - `:electorate` (default) — `Wᵢ/Nᵢ`, como no artigo, que calcula "a
    distribuição empírica de votos, Wᵢ/Nᵢ".
  - `:valid` — `Wᵢ/Vᵢ`, a convenção de boa parte das reimplementações. Muda a
    leitura do eixo y; use conscientemente.
"""
Base.@kwdef struct FingerprintConfig
    bins::Int = 50
    vote_axis::Symbol = :electorate
end

"""
    FingerprintResult

Resultado de [`election_fingerprint`](@ref). Implementa Tables.jl: uma linha por
casela ocupada do histograma, com `turnout`, `vote_rate` (centros da casela) e
`count`.

Campos: `turnout` e `vote_rate` por seção, `counts` (histograma 2-D),
`edges`, `m`, `n_excluded`, `nu` (taxa logarítmica de voto), `skewness`,
`kurtosis` e `vote_axis`.

`kurtosis` é a curtose **não** centrada em zero — a normal vale 3, como na
Fig. 3 do artigo.
"""
struct FingerprintResult
    turnout::Vector{Float64}
    vote_rate::Vector{Float64}
    counts::Matrix{Int}
    edges::LinRange{Float64,Int}
    m::Int
    n_excluded::Int
    nu::Vector{Float64}
    skewness::Float64
    kurtosis::Float64
    vote_axis::Symbol
end

"""
    log_vote_rate(winner_votes, electorate) -> Vector{Float64}

Taxa logarítmica de voto `νᵢ = log((Nᵢ − Wᵢ)/Wᵢ)` (Klimek et al. 2012, eq. da
seção "Data collapse"). Seções com `Wᵢ = 0` ou `Wᵢ ≥ Nᵢ` são omitidas — ν não
está definido nelas. O artigo chama essa exclusão de conservadora: descarta
justamente as seções de comparecimento e voto extremos, que são as mais
suspeitas.
"""
function log_vote_rate(winner_votes::AbstractVector{<:Integer},
                       electorate::AbstractVector{<:Integer})
    length(winner_votes) == length(electorate) ||
        throw(DimensionMismatch("winner_votes e electorate devem ter o mesmo comprimento"))
    ν = Float64[]
    sizehint!(ν, length(winner_votes))
    @inbounds for i in eachindex(winner_votes)
        W, N = Int(winner_votes[i]), Int(electorate[i])
        (W > 0 && W < N) || continue
        push!(ν, log((N - W) / W))
    end
    return ν
end

_skewness(v) = (μ = mean(v); s = std(v; corrected = false);
                s > 0 ? mean(((v .- μ) ./ s) .^ 3) : NaN)
_kurtosis(v) = (μ = mean(v); s = std(v; corrected = false);
                s > 0 ? mean(((v .- μ) ./ s) .^ 4) : NaN)

"""
    election_fingerprint(winner_votes, valid_votes, electorate; kwargs...)
        -> FingerprintResult

Constrói o fingerprint eleitoral de Klimek et al. (2012): o histograma 2-D do
número de seções por comparecimento `Vᵢ/Nᵢ` e por percentual de votos no
vencedor (`Wᵢ/Nᵢ` por default; veja [`FingerprintConfig`](@ref)).

Exige **três** entradas por seção — votos no vencedor, votos válidos e
eleitorado. É mais do que os outros testes do pacote pedem, e a maioria dos
repositórios de resultados não publica o eleitorado por seção.

Como ler:

- Um aglomerado único e compacto é o padrão de eleições sem alegação de fraude.
- Um borrão em direção ao canto superior direito indica manipulação
  incremental — cédulas retiradas de abstenções e da oposição e somadas ao
  vencedor.
- Um segundo aglomerado junto de (100 %, 100 %) indica manipulação extrema.

O resultado traz também `skewness` e `kurtosis` da taxa logarítmica de voto
([`log_vote_rate`](@ref)). Em eleições limpas elas ficam perto de (0, 3).

!!! warning "É um diagnóstico visual, não um teste"
    O fingerprint não produz p-valor, e assimetria/curtose isoladas não são
    teste de hipótese: em amostras grandes qualquer desvio ínfimo da
    normalidade rende uma estatística enorme. O artigo as usa como descrição
    comparativa entre países, e é assim que devem ser lidas. O modelo
    paramétrico de fraude `(fᵢ, f_e)` do artigo não está implementado — sua
    especificação completa está no Supporting Information, fora do preprint.
"""
function election_fingerprint(winner_votes::AbstractVector{<:Integer},
                              valid_votes::AbstractVector{<:Integer},
                              electorate::AbstractVector{<:Integer},
                              cfg::FingerprintConfig)
    m0 = length(winner_votes)
    (m0 == length(valid_votes) == length(electorate)) ||
        throw(DimensionMismatch("os três vetores devem ter o mesmo comprimento"))
    m0 > 0 || throw(ArgumentError("vetores vazios"))
    cfg.bins ≥ 2 || throw(ArgumentError("bins deve ser ≥ 2"))
    cfg.vote_axis in (:electorate, :valid) ||
        throw(ArgumentError("vote_axis deve ser :electorate ou :valid"))

    turnout = Float64[]; vote = Float64[]
    @inbounds for i in 1:m0
        W, V, N = Int(winner_votes[i]), Int(valid_votes[i]), Int(electorate[i])
        N > 0 || continue
        (0 ≤ W ≤ V ≤ N) || continue
        den = cfg.vote_axis === :electorate ? N : V
        den > 0 || continue
        push!(turnout, V / N)
        push!(vote, W / den)
    end
    m = length(turnout)
    m > 0 || throw(ArgumentError("nenhuma seção válida (exige 0 ≤ W ≤ V ≤ N, N > 0)"))

    edges = LinRange(0.0, 1.0, cfg.bins + 1)
    counts = zeros(Int, cfg.bins, cfg.bins)
    @inbounds for i in 1:m
        bx = clamp(floor(Int, turnout[i] * cfg.bins) + 1, 1, cfg.bins)
        by = clamp(floor(Int, vote[i] * cfg.bins) + 1, 1, cfg.bins)
        counts[bx, by] += 1
    end

    ν = log_vote_rate(winner_votes, electorate)
    sk = length(ν) > 2 ? _skewness(ν) : NaN
    ku = length(ν) > 3 ? _kurtosis(ν) : NaN

    return FingerprintResult(turnout, vote, counts, edges, m, m0 - m, ν, sk, ku,
                             cfg.vote_axis)
end

election_fingerprint(w::AbstractVector{<:Integer}, v::AbstractVector{<:Integer},
                     e::AbstractVector{<:Integer}; kwargs...) =
    election_fingerprint(w, v, e, FingerprintConfig(; kwargs...))

function Base.show(io::IO, ::MIME"text/plain", r::FingerprintResult)
    printstyled(io, "Fingerprint eleitoral — Klimek et al. (2012)\n"; bold = true)
    eixo = r.vote_axis === :electorate ? "W/N (eleitorado)" : "W/V (válidos)"
    println(io, "seções: $(r.m) (excluídas: $(r.n_excluded))   eixo y: $eixo")
    @printf(io, "comparecimento: mediana %.3f   votos no vencedor: mediana %.3f\n",
            median(r.turnout), median(r.vote_rate))
    @printf(io, "taxa log. de voto ν: n = %d   assimetria = %+.2f   curtose = %.2f\n",
            length(r.nu), r.skewness, r.kurtosis)
    println(io, "  (eleições limpas se agrupam em torno de assimetria 0, curtose 3)")

    # esboço do histograma em blocos, y crescendo para cima
    nb = size(r.counts, 1)
    passo = max(1, nb ÷ 20)
    mx = maximum(r.counts)
    if mx > 0
        rampa = (' ', '·', '░', '▒', '▓', '█')
        println(io, "  ┌", "─"^length(1:passo:nb), "┐ 100%")
        for by in reverse(1:passo:nb)
            print(io, "  │")
            for bx in 1:passo:nb
                s = sum(@view r.counts[bx:min(bx+passo-1, nb), by:min(by+passo-1, nb)])
                k = s == 0 ? 1 : clamp(ceil(Int, 5 * log1p(s) / log1p(mx)) + 1, 2, 6)
                print(io, rampa[k])
            end
            println(io, "│")
        end
        println(io, "  └", "─"^length(1:passo:nb), "┘ 0%")
        println(io, "   0%", " "^max(0, length(1:passo:nb) - 8), "comparecimento 100%")
    end
    # canto superior direito: assinatura de manipulação extrema
    q = max(1, ceil(Int, 0.95nb))
    extremo = sum(@view r.counts[q:nb, q:nb])
    if extremo > 0
        @printf(io, "seções no canto (≥95%% comparecimento e ≥95%% votos): %d (%.2f%%)\n",
                extremo, 100extremo / r.m)
    end
end
