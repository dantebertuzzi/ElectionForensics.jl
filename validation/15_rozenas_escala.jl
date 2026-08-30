# M5: a nula do teste de Rozenas superestima E[T] quando os shares se acumulam
# nas bordas. Diagnóstico: a banda global de Silverman (h≈0,07 sob Beta(.3,.3))
# é larga demais perto das fronteiras e empurra massa para dentro, onde as
# frações coarse são densas. Comparamos três remédios.
using ElectionForensics, Random, Distributions, Statistics, Printf
const FR = coarse_fractions(10)
const IDX = Dict{Rational{Int},Int}(f => j for (j, f) in enumerate(FR))

refl(p) = (while p < 0 || p > 1; p < 0 && (p = -p); p > 1 && (p = 2 - p); end; p)
logit(x) = log(x / (1 - x)); ilogit(z) = 1 / (1 + exp(-z))

"variantes do gerador de p̃ da nula"
function ptilde(kind, share, h, hl, pbar, lam, rng)
    if kind === :atual                       # Silverman + reflexão (patch atual)
        refl(share + h * randn(rng))
    elseif kind === :shrink                  # correção de variância (Silverman & Young 1987)
        refl(pbar + lam * (share - pbar) + h * randn(rng))
    elseif kind === :banda_menor             # banda fixa, só o bastante p/ borrar Farey
        refl(share + min(h, 0.02) * randn(rng))
    elseif kind === :logit                   # jitter na escala logit (respeita fronteiras)
        ilogit(logit(share) + hl * randn(rng))
    end
end

function teste(kind, votes, totals, B, rng)
    m = length(votes)
    T = count(i -> haskey(IDX, Int(votes[i]) // Int(totals[i])), 1:m)
    shares = Float64.(votes) ./ Float64.(totals)
    h, _ = ElectionForensics._silverman(shares)
    pbar = mean(shares); s2 = var(shares)
    v = mean(shares .* (1 .- shares) ./ totals)
    lam = sqrt(max(0.01, 1 - (h^2 + v) / max(s2, 1e-12)))
    # shares comprimidos p/ (0,1) abertos antes do logit
    sc = clamp.(shares, 1/(2*maximum(totals)), 1 - 1/(2*maximum(totals)))
    hl, _ = ElectionForensics._silverman(logit.(sc))
    hl = min(hl, 0.35)
    nulls = zeros(Int, B)
    for b in 1:B
        tb = 0
        for i in 1:m
            base = kind === :logit ? sc[i] : shares[i]
            p = ptilde(kind, base, h, hl, pbar, lam, rng)
            y = rand(rng, Binomial(Int(totals[i]), p))
            haskey(IDX, y // Int(totals[i])) && (tb += 1)
        end
        nulls[b] = tb
    end
    ((1 + count(≥(T), nulls)) / (B + 1), T, mean(nulls))
end

const R = 500
const VARIANTES = (:atual, :shrink, :banda_menor, :logit)

println("="^104)
println("M5 — ERRO TIPO I sob H0 (α = 0,05, m = 2000, B = 199).  Alvo: 0,050")
println("="^104)
@printf("%-14s", "dist. de p"); for k in VARIANTES; @printf("%-13s", k); end
@printf("%-11s%-11s\n", "E[T_obs]", "E[T*] atual")
for (nome, dist) in (("Beta(8,6)",Beta(8,6)), ("Beta(2,2)",Beta(2,2)),
                     ("Beta(1,1)",Beta(1,1)), ("Beta(.5,.5)",Beta(0.5,0.5)),
                     ("Beta(.3,.3)",Beta(0.3,0.3)), ("Beta(1.2,8)",Beta(1.2,8)))
    rej = Dict(k => zeros(R) for k in VARIANTES); to = zeros(R); tn = zeros(R)
    Threads.@threads for r in 1:R
        rng0 = Xoshiro(hash((nome, r)))
        totals = rand(rng0, 150:900, 2000)
        p = clamp.(rand(rng0, dist, 2000), 1e-9, 1-1e-9)
        votes = [rand(rng0, Binomial(totals[i], p[i])) for i in 1:2000]
        for k in VARIANTES
            pv, T, mn = teste(k, votes, totals, 199, Xoshiro(hash((nome, r, k))))
            rej[k][r] = pv
            k === :atual && (to[r] = T; tn[r] = mn)
        end
    end
    @printf("%-14s", nome); for k in VARIANTES; @printf("%-13.4f", mean(rej[k] .< 0.05)); end
    @printf("%-11.1f%-11.1f\n", mean(to), mean(tn))
end

println()
println("="^104)
println("M5 — PODER (metas redondas em ε das seções, shares Beta(2,2) polarizável)")
println("="^104)
@printf("%-10s%-10s", "dist.", "ε"); for k in VARIANTES; @printf("%-13s", k); end; println()
for (nome, dist) in (("Beta(8,6)",Beta(8,6)), ("Beta(.5,.5)",Beta(0.5,0.5)))
    for eps in (0.01, 0.02, 0.05)
        rej = Dict(k => zeros(R) for k in VARIANTES)
        Threads.@threads for r in 1:R
            rng0 = Xoshiro(hash((:pw, nome, eps, r)))
            totals = rand(rng0, 150:900, 2000)
            p = clamp.(rand(rng0, dist, 2000), 1e-9, 1-1e-9)
            votes = [rand(rng0, Binomial(totals[i], p[i])) for i in 1:2000]
            alvos = (1//2, 3//5, 2//3, 7//10, 3//4, 4//5)
            for i in 1:round(Int, eps*2000)
                f = rand(rng0, alvos); votes[i] = (totals[i]*numerator(f))÷denominator(f)
            end
            for k in VARIANTES
                rej[k][r] = teste(k, votes, totals, 199, Xoshiro(hash((:pw,nome,eps,r,k))))[1]
            end
        end
        @printf("%-10s%-10.2f", nome, eps); for k in VARIANTES; @printf("%-13.3f", mean(rej[k].<0.05)); end; println()
    end
end
