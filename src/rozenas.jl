# ── Rozenas (2017): excesso de frações coarse ────────────────────────

"""
    coarse_fractions(max_denom::Int = 10) -> Vector{Rational{Int}}

Frações irredutíveis k/d com 2 ≤ d ≤ `max_denom` e 0 < k/d < 1
(sequência de Farey sem os extremos). São os percentuais "redondos"
(1/2, 2/3, 3/4, ...) nos quais resultados fabricados tendem a se
concentrar.
"""
function coarse_fractions(max_denom::Int = 10)
    max_denom ≥ 2 || throw(ArgumentError("max_denom deve ser ≥ 2"))
    fr = [k // d for d in 2:max_denom for k in 1:(d - 1) if gcd(k, d) == 1]
    return sort!(fr)
end

struct RozenasResult
    fractions::Vector{Rational{Int}}
    observed::Vector{Int}
    null_mean::Vector{Float64}
    null_sd::Vector{Float64}
    total_observed::Int
    null_totals::Vector{Int}
    pvalue::Float64
    zscore::Float64
    m::Int
    B::Int
    h::Float64
end

"""
    rozenas_test(votes::AbstractVector{<:Integer},
                 totals::AbstractVector{<:Integer};
                 fractions = coarse_fractions(10),
                 B::Int = 999,
                 h::Union{Nothing,Real} = nothing,
                 rng::AbstractRNG = Random.default_rng())
        -> RozenasResult

Teste de frações coarse de Rozenas (2017, *Political Analysis*): conta
quantas seções têm percentual de votos `votes[i] / totals[i]` **exatamente**
igual a uma fração redonda k/d (comparação em aritmética racional, sem
erro de ponto flutuante) e compara com a distribuição nula obtida por
bootstrap paramétrico.

Distribuição nula: em cada uma das `B` réplicas, o percentual de cada
seção é perturbado por um kernel gaussiano de banda `h` (default:
Silverman sobre os percentuais observados — aproximação da densidade
kernel reamostrada do artigo) e novas contagens são sorteadas de
`Binomial(totals[i], p̃ᵢ)`. Isso estima quantas seções cairiam em frações
coarse *por acaso*, dado o tamanho das seções e o formato da distribuição
de percentuais.

Estatística de teste: `T = Σ` seções em frações coarse. O p-valor é
unilateral, `P(T* ≥ T_obs)`, com correção `(1 + #{T*_b ≥ T}) / (B + 1)`.

Notas:
- Percentuais 0 e 1 (seções unânimes) não contam como frações coarse.
- Perturbar em torno do percentual observado torna o teste levemente
  conservador sob fraude maciça (a nula "herda" parte do excesso);
  um resultado significativo é, portanto, evidência forte.
"""
function rozenas_test(votes::AbstractVector{<:Integer},
                      totals::AbstractVector{<:Integer};
                      fractions::Vector{Rational{Int}} = coarse_fractions(10),
                      B::Int = 999,
                      h::Union{Nothing,Real} = nothing,
                      rng::AbstractRNG = Random.default_rng())
    m = length(votes)
    m == length(totals) ||
        throw(DimensionMismatch("votes e totals devem ter o mesmo comprimento"))
    m > 0 || throw(ArgumentError("vetores vazios"))
    B ≥ 99 || throw(ArgumentError("B deve ser ≥ 99"))
    @inbounds for i in 1:m
        totals[i] > 0 || throw(ArgumentError("totals[$i] deve ser > 0"))
        0 ≤ votes[i] ≤ totals[i] ||
            throw(ArgumentError("votes[$i] fora de [0, totals[$i]]"))
    end

    J = length(fractions)
    frac_index = Dict{Rational{Int},Int}(f => j for (j, f) in enumerate(fractions))

    # contagem observada (aritmética exata via Rational)
    observed = zeros(Int, J)
    @inbounds for i in 1:m
        r = Int(votes[i]) // Int(totals[i])
        j = get(frac_index, r, 0)
        j > 0 && (observed[j] += 1)
    end
    T_obs = sum(observed)

    # nula por bootstrap paramétrico com jitter gaussiano
    shares = Float64.(votes) ./ Float64.(totals)
    hval = h === nothing ? _silverman(shares) : Float64(h)
    hval > 0 || throw(ArgumentError("h deve ser > 0"))

    null_counts = zeros(Int, B, J)
    null_totals = zeros(Int, B)
    @inbounds for b in 1:B
        tb = 0
        for i in 1:m
            p = clamp(shares[i] + hval * randn(rng), 0.0, 1.0)
            y = rand(rng, Binomial(Int(totals[i]), p))
            r = y // Int(totals[i])
            j = get(frac_index, r, 0)
            if j > 0
                null_counts[b, j] += 1
                tb += 1
            end
        end
        null_totals[b] = tb
    end

    null_mean = vec(mean(null_counts; dims = 1))
    null_sd = vec(std(null_counts; dims = 1))
    μT = mean(null_totals)
    σT = std(null_totals)
    zscore = σT > 0 ? (T_obs - μT) / σT : (T_obs > μT ? Inf : 0.0)
    pvalue = (1 + count(≥(T_obs), null_totals)) / (B + 1)

    return RozenasResult(fractions, observed, null_mean, null_sd,
                         T_obs, null_totals, pvalue, zscore, m, B, hval)
end

function Base.show(io::IO, ::MIME"text/plain", r::RozenasResult)
    printstyled(io, "Teste de frações coarse — Rozenas (2017)\n"; bold = true)
    @printf(io, "seções: %d   réplicas: %d   h (jitter): %.4f\n", r.m, r.B, r.h)
    @printf(io, "T_obs = %d   E[T|H₀] = %.1f   z = %.2f   p = %s\n",
            r.total_observed, mean(r.null_totals), r.zscore, _fmt_p(r.pvalue))

    # frações mais anômalas (excesso padronizado)
    zs = [r.null_sd[j] > 0 ? (r.observed[j] - r.null_mean[j]) / r.null_sd[j] :
          (r.observed[j] > r.null_mean[j] ? Inf : 0.0) for j in eachindex(r.fractions)]
    ord = sortperm(zs; rev = true)
    top = [j for j in ord if r.observed[j] > 0][1:min(10, count(>(0), r.observed))]
    if !isempty(top)
        println(io, "frações com maior excesso:")
        println(io, "  fração    obs   esp     z")
        for j in top
            @printf(io, "  %-8s %4d  %6.1f  %5.2f\n",
                    string(numerator(r.fractions[j]), "/", denominator(r.fractions[j])),
                    r.observed[j], r.null_mean[j], zs[j])
        end
    end
    print(io, "H₀ (sem excesso em frações coarse): ")
    printstyled(io, r.pvalue < 0.05 ? "rejeitada" : "não rejeitada";
                color = _sig_color(r.pvalue), bold = true)
    println(io)
end
