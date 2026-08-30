# ── mistura de Beta-Binomiais (nula de Rozenas 2017) ─────────────────
#
# Reimplementação do modelo do pacote `spikes` do próprio Rozenas (CRAN, 2016).
# Em vez de perturbar os percentuais observados com um kernel — que herda o
# ruído binomial e precisa de tratamento ad hoc nas fronteiras — ajusta-se uma
# mistura de Beta-Binomiais à distribuição LATENTE dos percentuais e reamostra-se
# da posterior de cada seção:
#
#     p̃ᵢ ~ Beta(yᵢ + αₖ, nᵢ − yᵢ + βₖ),   k sorteado ∝ responsabilidade
#     yᵢ* ~ Binomial(nᵢ, p̃ᵢ)
#
# A posterior deconvolui o ruído binomial (seções pequenas encolhem para a
# prior) e vive naturalmente em [0,1], então não há fronteira para tratar.
#
# Duas divergências deliberadas em relação ao `spikes`, ambas defeitos que a
# revisão identificou no original:
#
#  1. Seleção do número de componentes por BIC, não pela correlação entre
#     densidades de L consecutivos. O critério original não tem interpretação
#     estatística e depende da ordem em que L cresce.
#  2. Todos os laços têm teto de iterações. O `spikes` usa `while (delta < 0.99)`
#     sem limite, que pode não terminar.

"""
    BetaBinomialMixture

Mistura de Beta-Binomiais ajustada por EM, usada como nula do teste de Rozenas
quando `null = :betabinomial`.

Campos: `weights`, `alpha`, `beta` (um por componente), `loglik` (núcleo da
log-verossimilhança — omite o termo `log C(n,y)`, constante em α, β e L),
`bic`, `components`, `iterations`, e o diagnóstico de aderência `misfit` /
`misfit_range`.

`misfit` é a maior razão observado/predito entre as faixas de percentual com
massa não desprezível, e `misfit_range` diz em qual faixa ela ocorre. Sendo
paramétrica, a mistura pode ajustar bem o corpo da distribuição e errar numa
região fina — e é exatamente ali que a nula passa a inventar excesso. Nos dados
russos de 2012 a cauda esquerda é subestimada em 3,8×, o que fabrica dez
frações "significativas" entre 1/7 e 1/4 (`validation/21_nula_betabinomial.jl`).
Um `misfit` acima de ~2 é motivo para preferir `null = :kernel`.
"""
struct BetaBinomialMixture
    weights::Vector{Float64}
    alpha::Vector{Float64}
    beta::Vector{Float64}
    loglik::Float64
    bic::Float64
    components::Int
    iterations::Int
    misfit::Float64
    misfit_range::Tuple{Float64,Float64}
end

# Núcleo do log-pmf Beta-Binomial. O termo log C(n,y) é constante em (α,β) e em
# L, então cancela nas responsabilidades, no passo M e na comparação por BIC.
@inline function _bb_logkernel(y::Int, n::Int, a::Float64, b::Float64)
    return logbeta(y + a, n - y + b) - logbeta(a, b)
end

"""
Passo E: preenche `W` com as responsabilidades e devolve o núcleo da
log-verossimilhança. Feito em escala log com `logsumexp` — com milhares de
seções, o produto direto das densidades sofre underflow.
"""
function _bb_estep!(W::Matrix{Float64}, y::Vector{Int}, n::Vector{Int},
                    weights, alpha, beta)
    m, L = length(y), length(weights)
    logw = log.(weights)
    ll = 0.0
    @inbounds for i in 1:m
        mx = -Inf
        for k in 1:L
            v = logw[k] + _bb_logkernel(y[i], n[i], alpha[k], beta[k])
            W[i, k] = v
            v > mx && (mx = v)
        end
        if !isfinite(mx)                     # componente degenerado
            for k in 1:L; W[i, k] = 1 / L; end
            continue
        end
        s = 0.0
        for k in 1:L; s += exp(W[i, k] - mx); end
        lse = mx + log(s)
        ll += lse
        for k in 1:L; W[i, k] = exp(W[i, k] - lse); end
    end
    return ll
end

# Nelder–Mead em 2-D com caixa, para o passo M. Duas dimensões e superfície
# suave: não compensa arrastar uma dependência de otimização.
function _nelder_mead(f, x0::NTuple{2,Float64}, lo::Float64, hi::Float64;
                      maxiter::Int = 200, tol::Float64 = 1e-9)
    clamp2(p) = (clamp(p[1], lo, hi), clamp(p[2], lo, hi))
    step = 0.5 * (hi - lo) / 10
    simplex = [clamp2(x0),
               clamp2((x0[1] + step, x0[2])),
               clamp2((x0[1], x0[2] + step))]
    fv = [f(p) for p in simplex]
    for _ in 1:maxiter
        ord = sortperm(fv)
        simplex, fv = simplex[ord], fv[ord]
        (abs(fv[end] - fv[1]) ≤ tol * (abs(fv[1]) + tol)) && break
        c = ((simplex[1][1] + simplex[2][1]) / 2, (simplex[1][2] + simplex[2][2]) / 2)
        w = simplex[3]
        refl = clamp2((c[1] + (c[1] - w[1]), c[2] + (c[2] - w[2])))
        fr = f(refl)
        if fr < fv[1]
            exp_ = clamp2((c[1] + 2(c[1] - w[1]), c[2] + 2(c[2] - w[2])))
            fe = f(exp_)
            simplex[3], fv[3] = fe < fr ? (exp_, fe) : (refl, fr)
        elseif fr < fv[2]
            simplex[3], fv[3] = refl, fr
        else
            con = clamp2((c[1] + 0.5(w[1] - c[1]), c[2] + 0.5(w[2] - c[2])))
            fc = f(con)
            if fc < fv[3]
                simplex[3], fv[3] = con, fc
            else
                for j in 2:3
                    simplex[j] = clamp2((simplex[1][1] + 0.5(simplex[j][1] - simplex[1][1]),
                                         simplex[1][2] + 0.5(simplex[j][2] - simplex[1][2])))
                    fv[j] = f(simplex[j])
                end
            end
        end
    end
    i = argmin(fv)
    return simplex[i], fv[i]
end

"Ajusta uma mistura de `L` Beta-Binomiais por EM. Devolve `nothing` se degenerar."
function _fit_bb_mixture(y::Vector{Int}, n::Vector{Int}, L::Int;
                         maxiter::Int, tol::Float64, lo::Float64, hi::Float64,
                         rng::AbstractRNG)
    m = length(y)
    shares = [y[i] / n[i] for i in 1:m]
    # Inicialização por momentos: quantis dos percentuais observados dão a
    # média de cada componente; a concentração vem da dispersão amostral.
    μ0 = L == 1 ? [clamp(mean(shares), 0.02, 0.98)] :
         [clamp(quantile(shares, (2k - 1) / (2L)), 0.02, 0.98) for k in 1:L]
    s2 = max(var(shares), 1e-4)
    conc = max(mean(μ0) * (1 - mean(μ0)) / s2 - 1, 1.0)
    alpha = [clamp(μ0[k] * conc, exp(lo), exp(hi)) for k in 1:L]
    beta  = [clamp((1 - μ0[k]) * conc, exp(lo), exp(hi)) for k in 1:L]
    weights = fill(1 / L, L)
    # perturbação para quebrar simetria entre componentes iguais
    @inbounds for k in 1:L
        alpha[k] = clamp(alpha[k] * exp(0.1 * randn(rng)), exp(lo), exp(hi))
        beta[k]  = clamp(beta[k]  * exp(0.1 * randn(rng)), exp(lo), exp(hi))
    end

    W = Matrix{Float64}(undef, m, L)
    ll_old = -Inf; ll = -Inf; iters = 0
    for it in 1:maxiter
        iters = it
        ll = _bb_estep!(W, y, n, weights, alpha, beta)
        isfinite(ll) || return nothing
        # passo M
        @inbounds for k in 1:L
            wk = @view W[:, k]
            sw = sum(wk)
            weights[k] = sw / m
            sw < 1e-8 && continue            # componente vazio: mantém como está
            # maximiza a verossimilhança ponderada em (log α, log β)
            obj = function (p)
                a, b = exp(p[1]), exp(p[2])
                s = 0.0
                for i in 1:m
                    wi = wk[i]
                    wi > 1e-12 && (s += wi * _bb_logkernel(y[i], n[i], a, b))
                end
                return -s                     # Nelder–Mead minimiza
            end
            # EM generalizado: basta melhorar a verossimilhança a cada passo,
            # não maximizá-la — daí o teto baixo de iterações do Nelder–Mead.
            best, _ = _nelder_mead(obj, (log(alpha[k]), log(beta[k])), lo, hi;
                                   maxiter = 60)
            alpha[k], beta[k] = exp(best[1]), exp(best[2])
        end
        sw = sum(weights)
        sw > 0 ? (weights ./= sw) : return nothing
        # convergência só a partir da segunda iteração: com `ll_old = -Inf` a
        # tolerância relativa é infinita e o laço pararia de imediato.
        if it > 1 && abs(ll - ll_old) ≤ tol * (abs(ll_old) + tol)
            break
        end
        ll_old = ll
    end
    npar = 3L - 1
    mis, faixa = _bb_goodness_of_fit(y, n, weights, alpha, beta, W)
    return BetaBinomialMixture(copy(weights), copy(alpha), copy(beta),
                               ll, -2ll + npar * log(m), L, iters, mis, faixa)
end

"""
Checagem preditiva a posteriori: gera uma réplica da nula e compara a
distribuição de percentuais com a observada, faixa a faixa. Devolve a maior
razão observado/predito entre as faixas com massa ≥ `piso`, e a faixa onde
ocorre. Sem isso, um desajuste numa região fina passa despercebido e vira
"excesso" nas frações que caem ali.
"""
function _bb_goodness_of_fit(y, n, weights, alpha, beta, W; nbins::Int = 20,
                             piso::Float64 = 5e-4)
    m = length(y)
    mix = BetaBinomialMixture(collect(weights), collect(alpha), collect(beta),
                              NaN, NaN, length(weights), 0, NaN, (NaN, NaN))
    buf = Vector{Int}(undef, m)
    _resample_betabinomial!(buf, y, n, mix, W, MersenneTwister(0x5eed))
    obs = zeros(Int, nbins); nul = zeros(Int, nbins)
    bin(x) = clamp(floor(Int, x * nbins) + 1, 1, nbins)
    @inbounds for i in 1:m
        obs[bin(y[i] / n[i])] += 1
        nul[bin(buf[i] / n[i])] += 1
    end
    pior = 1.0; faixa = (0.0, 1.0)
    for b in 1:nbins
        po = obs[b] / m
        po < piso && continue
        pn = max(nul[b] / m, 1 / m)
        razao = max(po / pn, pn / po)
        if razao > pior
            pior = razao
            faixa = ((b - 1) / nbins, b / nbins)
        end
    end
    return pior, faixa
end

"""
    fit_betabinomial_mixture(votes, totals; max_components = 5, kwargs...)
        -> BetaBinomialMixture

Ajusta por EM a mistura de Beta-Binomiais que descreve a distribuição latente
de percentuais, escolhendo o número de componentes por **BIC** entre 1 e
`max_components`. É a nula usada por `rozenas_test(...; null = :betabinomial)`.

Segue o modelo do pacote `spikes` de Rozenas (CRAN, 2016), com duas
divergências deliberadas: seleção por BIC em vez da correlação entre densidades
de L consecutivos, e teto de iterações em todos os laços.
"""
function fit_betabinomial_mixture(votes::AbstractVector{<:Integer},
                                  totals::AbstractVector{<:Integer};
                                  max_components::Int = 5,
                                  maxiter::Int = 100,
                                  tol::Float64 = 1e-7,
                                  log_bounds::Tuple{Real,Real} = (-3.0, 9.0),
                                  rng::AbstractRNG = Random.default_rng())
    max_components ≥ 1 || throw(ArgumentError("max_components deve ser ≥ 1"))
    lo, hi = Float64(log_bounds[1]), Float64(log_bounds[2])
    lo < hi || throw(ArgumentError("log_bounds inválido"))
    y = [Int(v) for v in votes]; n = [Int(t) for t in totals]
    melhor = nothing
    for L in 1:max_components
        cand = _fit_bb_mixture(y, n, L; maxiter = maxiter, tol = tol,
                               lo = lo, hi = hi, rng = rng)
        cand === nothing && continue
        (melhor === nothing || cand.bic < melhor.bic) && (melhor = cand)
    end
    melhor === nothing &&
        throw(ArgumentError("não foi possível ajustar a mistura Beta-Binomial"))
    return melhor
end

"""
Reamostra as contagens da nula: para cada seção sorteia um componente ∝
responsabilidade, tira o percentual latente da posterior `Beta(y+α, n−y+β)` e
redesenha a contagem de `Binomial(n, p̃)`.
"""
function _resample_betabinomial!(out::Vector{Int}, y::Vector{Int}, n::Vector{Int},
                                 mix::BetaBinomialMixture, W::Matrix{Float64},
                                 rng::AbstractRNG)
    L = mix.components
    @inbounds for i in eachindex(y)
        k = 1
        if L > 1
            u = rand(rng); acc = 0.0
            for j in 1:L
                acc += W[i, j]
                if u ≤ acc; k = j; break; end
                k = j
            end
        end
        p = rand(rng, Beta(y[i] + mix.alpha[k], n[i] - y[i] + mix.beta[k]))
        out[i] = rand(rng, Binomial(n[i], p))
    end
    return out
end
