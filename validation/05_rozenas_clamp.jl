# O clamp(share + h·randn, 0, 1) empilha massa EXATAMENTE em 0 e 1.
# Nesses pontos Binomial(n, 0) = 0 e Binomial(n, 1) = n, cujas frações (0 e 1)
# nunca são coarse -> a nula subconta T* -> p-valor deflacionado.
# Comparamos com reflexão na fronteira (que preserva a densidade).
using ElectionForensics, Random, Distributions, Statistics, Printf

# ── reimplementação idêntica ao pacote, exceto o tratamento da fronteira ──
function rozenas_T(votes, totals, fractions, B, h, rng; boundary = :clamp)
    idx = Dict{Rational{Int},Int}(f => j for (j, f) in enumerate(fractions))
    T_obs = count(i -> haskey(idx, Int(votes[i]) // Int(totals[i])), eachindex(votes))
    shares = Float64.(votes) ./ Float64.(totals)
    nulls = zeros(Int, B)
    for b in 1:B
        tb = 0
        for i in eachindex(shares)
            p = shares[i] + h * randn(rng)
            if boundary === :clamp
                p = clamp(p, 0.0, 1.0)
            else                                   # reflexão em 0 e 1
                while p < 0 || p > 1
                    p < 0 && (p = -p)
                    p > 1 && (p = 2 - p)
                end
            end
            y = rand(rng, Binomial(Int(totals[i]), p))
            haskey(idx, y // Int(totals[i])) && (tb += 1)
        end
        nulls[b] = tb
    end
    ((1 + count(≥(T_obs), nulls)) / (B + 1), T_obs, mean(nulls))
end

const R = 800
const FR = coarse_fractions(10)

println("="^96)
println("ERRO TIPO I sob H0 em eleições POLARIZADAS (muitas seções quase unânimes)")
println("m = 2000, B = 199, α = 0.05.  clamp = comportamento atual do pacote")
println("="^96)
@printf("%-16s %-9s %-9s %-11s %-11s %-9s %-9s\n",
        "dist. de p", "rej clamp", "rej reflex", "E[T*] clamp", "E[T*] reflex", "E[T_obs]", "% clamp")

for (nome, dist) in (("Beta(2,2)",   Beta(2,2)),
                     ("Beta(1,1)",   Beta(1,1)),
                     ("Beta(.5,.5)", Beta(0.5,0.5)),
                     ("Beta(.3,.3)", Beta(0.3,0.3)),
                     ("Beta(1.2,8)", Beta(1.2,8)))
    pc = zeros(R); pr = zeros(R); tc = zeros(R); tr_ = zeros(R); to = zeros(R)
    Threads.@threads for r in 1:R
        rng = Xoshiro(hash((nome, r)))
        totals = rand(rng, 150:900, 2000)
        p = clamp.(rand(rng, dist, 2000), 1e-9, 1-1e-9)
        votes = [rand(rng, Binomial(totals[i], p[i])) for i in 1:2000]
        h = ElectionForensics._silverman(votes ./ totals)
        a = rozenas_T(votes, totals, FR, 199, h, Xoshiro(hash((nome, r, :c))); boundary = :clamp)
        b = rozenas_T(votes, totals, FR, 199, h, Xoshiro(hash((nome, r, :r))); boundary = :reflect)
        pc[r], to[r], tc[r] = a
        pr[r], _,     tr_[r] = b
    end
    # fração clampada
    rng = Xoshiro(hash((nome, 1)))
    totals = rand(rng, 150:900, 2000); p = clamp.(rand(rng, dist, 2000), 1e-9, 1-1e-9)
    votes = [rand(rng, Binomial(totals[i], p[i])) for i in 1:2000]
    h = ElectionForensics._silverman(votes ./ totals)
    nc = mean([(q = votes[i]/totals[i] + h*randn(rng); q <= 0 || q >= 1) for _ in 1:200 for i in 1:2000])
    @printf("%-16s %-9.4f %-9.4f %-11.1f %-11.1f %-9.1f %-9.2f\n",
            nome, mean(pc .< 0.05), mean(pr .< 0.05), mean(tc), mean(tr_), mean(to), 100nc)
end
