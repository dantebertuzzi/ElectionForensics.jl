# ── Rozenas (2017): excesso de frações coarse ────────────────────────

# Fração y/n reduzida como par de inteiros — ver `frac_index` em `rozenas_test`.
@inline function _reduced(y::Int, n::Int)
    g = gcd(y, n)
    return g == 0 ? (0, 0) : (y ÷ g, n ÷ g)
end

_logit(x) = log(x / (1 - x))
_ilogit(z) = 1 / (1 + exp(-z))

# Correção de continuidade por seção: leva shares de 0 e 1 para dentro de (0,1)
# sem colapsar seções distintas no mesmo ponto.
_compress(y, n) = (y + 0.5) / (n + 1)

# Duas escalas para o jitter da nula:
#
# `:logit` (default) — jitter gaussiano em log(p/(1-p)). Perto de 0 e 1 o passo
#   em escala de probabilidade encolhe sozinho, então a nula não empurra massa
#   das seções quase unânimes para o miolo, onde as frações coarse são densas.
#   Sob shares em U isso levava o erro tipo I a 0,000 (α = 0,05) e custava
#   poder: 0,48 contra 0,67 a ε = 2 % (validation/15_rozenas_escala.jl).
#
# `:reflect` — jitter na escala de probabilidade com reflexão nas fronteiras.
#   Preserva o comportamento histórico; `h` passa a ser lido nessa escala.
function _jitter(boundary::Symbol, i::Int, shares, jshares, h::Float64, rng)
    if boundary === :logit
        return _ilogit(_logit(jshares[i]) + h * randn(rng))
    else
        p = shares[i] + h * randn(rng)
        while p < 0.0 || p > 1.0
            p < 0.0 && (p = -p)
            p > 1.0 && (p = 2.0 - p)
        end
        return p
    end
end

# Teto da banda em escala logit: sem ele, uma distribuição de shares muito
# dispersa produziria jitter que apaga o formato junto com a estrutura de Farey.
const _HLOGIT_CAP = 0.35

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

"""
    RozenasResult

Resultado de [`rozenas_test`](@ref).

Campos: `fractions` (as frações coarse testadas), `observed` (seções em cada
fração), `null_mean`/`null_sd` (nula por fração), `total_observed` (estatística
agregada `T`), `null_totals` (as `B` réplicas de `T`), `pvalue` (Monte Carlo,
unilateral), `zscore`, `qvalues` (p-valores por fração já ajustados por
Benjamini–Hochberg), `m` (seções), `B` (réplicas) e `h` (banda do jitter).

Só `pvalue` é um teste; os excessos por fração são exploratórios e por isso
vêm com `qvalues` — sob H₀, o maior z entre as 31 frações tem mediana 2,4
(`validation/11_zscores_selecao.jl`).
"""
struct RozenasResult
    fractions::Vector{Rational{Int}}
    observed::Vector{Int}
    null_mean::Vector{Float64}
    null_sd::Vector{Float64}
    total_observed::Int
    null_totals::Vector{Int}
    pvalue::Float64
    zscore::Float64
    qvalues::Vector{Float64}
    m::Int
    B::Int
    h::Float64
    boundary::Symbol
    null::Symbol
    mixture::Union{Nothing,BetaBinomialMixture}
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
                      boundary::Symbol = :logit,
                      null::Symbol = :kernel,
                      max_components::Int = 5,
                      rng::AbstractRNG = Random.default_rng())
    boundary in (:logit, :reflect) ||
        throw(ArgumentError("boundary deve ser :logit ou :reflect"))
    null in (:kernel, :betabinomial) ||
        throw(ArgumentError("null deve ser :kernel ou :betabinomial"))
    m = length(votes)
    m == length(totals) ||
        throw(DimensionMismatch("votes e totals devem ter o mesmo comprimento"))
    m > 0 || throw(ArgumentError("vetores vazios"))
    m ≥ 30 || @warn "poucas seções para o teste de Rozenas (m = $m); \
        o p-valor Monte Carlo tem resolução muito baixa"
    B ≥ 99 || throw(ArgumentError("B deve ser ≥ 99"))
    @inbounds for i in 1:m
        totals[i] > 0 || throw(ArgumentError("totals[$i] deve ser > 0"))
        0 ≤ votes[i] ≤ totals[i] ||
            throw(ArgumentError("votes[$i] fora de [0, totals[$i]]"))
    end

    J = length(fractions)
    # Chave (num, den) reduzida em vez de `Rational`: `hash(::Rational)` passa
    # por `ldexp`/`reinterpret` para casar com o hash de Float, o que custa caro
    # nas m×B buscas do laço da nula. A semântica é idêntica.
    frac_index = Dict{Tuple{Int,Int},Int}(
        (numerator(f), denominator(f)) => j for (j, f) in enumerate(fractions))

    # contagem observada (aritmética exata via Rational)
    observed = zeros(Int, J)
    @inbounds for i in 1:m
        j = get(frac_index, _reduced(Int(votes[i]), Int(totals[i])), 0)
        j > 0 && (observed[j] += 1)
    end
    T_obs = sum(observed)

    # nula por bootstrap paramétrico com jitter gaussiano
    shares = Float64.(votes) ./ Float64.(totals)
    # shares comprimidos para (0,1) aberto — só o gerador da nula os usa;
    # a estatística observada continua vindo da aritmética racional exata.
    jshares = [_compress(Float64(votes[i]), Float64(totals[i])) for i in 1:m]

    if h === nothing
        hval, degenerado = boundary === :logit ? _silverman(_logit.(jshares)) :
                                                 _silverman(shares)
        boundary === :logit && (hval = min(hval, _HLOGIT_CAP))
        degenerado && @warn "os percentuais observados não têm dispersão; \
            usando banda mínima h = $(_H_FLOOR). A nula é pouco informativa."
    else
        hval = Float64(h)
    end
    hval > 0 || throw(ArgumentError("h deve ser > 0"))

    # Ajuste da mistura Beta-Binomial, quando pedida. Feito uma única vez: as
    # B réplicas reamostram da mesma posterior.
    mixture = nothing
    Wresp = zeros(Float64, 0, 0)
    yv = Int[]; nv = Int[]
    if null === :betabinomial
        yv = [Int(votes[i]) for i in 1:m]
        nv = [Int(totals[i]) for i in 1:m]
        mixture = fit_betabinomial_mixture(yv, nv;
                                           max_components = max_components, rng = rng)
        if mixture.misfit > 2.0
            @warn "a mistura Beta-Binomial não descreve bem a faixa de \
                percentuais [$(round(mixture.misfit_range[1], digits=2)), \
                $(round(mixture.misfit_range[2], digits=2))): observado é \
                $(round(mixture.misfit, digits=1))× o predito. A nula vai \
                inventar excesso nas frações dessa faixa — prefira null = :kernel."
        end
        Wresp = Matrix{Float64}(undef, m, mixture.components)
        _bb_estep!(Wresp, yv, nv, mixture.weights, mixture.alpha, mixture.beta)
    end

    null_counts = zeros(Int, B, J)
    null_totals = zeros(Int, B)
    buf = Vector{Int}(undef, m)
    @inbounds for b in 1:B
        if null === :betabinomial
            _resample_betabinomial!(buf, yv, nv, mixture, Wresp, rng)
        else
            for i in 1:m
                p = _jitter(boundary, i, shares, jshares, hval, rng)
                buf[i] = rand(rng, Binomial(Int(totals[i]), p))
            end
        end
        tb = 0
        for i in 1:m
            j = get(frac_index, _reduced(buf[i], Int(totals[i])), 0)
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

    # p-valores por fração (Monte Carlo, unilateral) ajustados por
    # Benjamini–Hochberg: sem isso, o maior z entre 31 frações passa de 2 em
    # 68 % das eleições limpas (validation/11_zscores_selecao.jl).
    praw = [(1 + count(≥(observed[j]), @view null_counts[:, j])) / (B + 1) for j in 1:J]
    qvalues = _bh_adjust(praw)

    return RozenasResult(fractions, observed, null_mean, null_sd,
                         T_obs, null_totals, pvalue, zscore, qvalues, m, B,
                         hval, boundary, null, mixture)
end

function Base.show(io::IO, ::MIME"text/plain", r::RozenasResult)
    printstyled(io, "Teste de frações coarse — Rozenas (2017)\n"; bold = true)
    if r.null === :betabinomial
        @printf(io, "seções: %d   réplicas: %d   nula: mistura Beta-Binomial (%d componentes)\n",
                r.m, r.B, r.mixture.components)
    else
        @printf(io, "seções: %d   réplicas: %d   nula: kernel (escala %s, h = %.4f)\n",
                r.m, r.B, r.boundary === :logit ? "logit" : "prob.", r.h)
    end
    @printf(io, "T_obs = %d   E[T|H₀] = %.1f   z = %.2f   p = %s\n",
            r.total_observed, mean(r.null_totals), r.zscore, _fmt_p(r.pvalue))

    # frações mais anômalas — exploratório, com q-valor de Benjamini–Hochberg
    zs = [r.null_sd[j] > 0 ? (r.observed[j] - r.null_mean[j]) / r.null_sd[j] :
          (r.observed[j] > r.null_mean[j] ? Inf : 0.0) for j in eachindex(r.fractions)]
    ord = sortperm(zs; rev = true)
    top = [j for j in ord if r.observed[j] > 0][1:min(10, count(>(0), r.observed))]
    if !isempty(top)
        println(io, "frações com maior excesso (exploratório; q = BH sobre $(length(r.fractions)) frações):")
        println(io, "  fração    obs   esp     z      q")
        for j in top
            @printf(io, "  %-8s %4d  %6.1f  %5.2f  %s\n",
                    string(numerator(r.fractions[j]), "/", denominator(r.fractions[j])),
                    r.observed[j], r.null_mean[j], zs[j], _fmt_p(r.qvalues[j]))
        end
    end
    print(io, "H₀ (sem excesso em frações coarse): ")
    printstyled(io, r.pvalue < 0.05 ? "rejeitada" : "não rejeitada";
                color = _sig_color(r.pvalue), bold = true)
    println(io)
end

function rozenas_test(votes::AbstractVector, totals::AbstractVector;
                      skipmissing::Bool = false, kwargs...)
    (votes isa AbstractVector{<:Integer} && totals isa AbstractVector{<:Integer}) &&
        throw(MethodError(rozenas_test, (votes, totals)))
    skipmissing && (any(ismissing, votes) || any(ismissing, totals)) &&
        throw(ArgumentError("`skipmissing` descartaria seções isoladamente e \
            desalinharia votes/totals; remova as seções incompletas antes"))
    return rozenas_test(_coerce_counts(votes; name = "votes"),
                        _coerce_counts(totals; name = "totals"); kwargs...)
end
