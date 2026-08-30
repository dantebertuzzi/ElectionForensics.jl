# Calibração sob H0 do teste de frações coarse.
# DGP limpo EXATAMENTE no modelo que o teste pressupõe:
#   p_i ~ F (suave),  votes_i ~ Binomial(totals_i, p_i)
using ElectionForensics, Random, Distributions, Statistics, HypothesisTests, Printf

const R = 1000
const B = 199

function replica(rng, m, totrange, dist)
    totals = rand(rng, totrange, m)
    p = clamp.(rand(rng, dist, m), 0.001, 0.999)
    votes = [rand(rng, Binomial(totals[i], p[i])) for i in 1:m]
    return votes, totals
end

println("="^92)
println("ROZENAS — erro tipo I sob H0 (R = $R réplicas, B = $B)")
println("="^92)
@printf("%-22s %-7s %-8s %-8s %-8s %-8s %-8s %-8s\n",
        "cenário", "m", "KS p", "rej@.01", "rej@.05", "rej@.10", "E[T_obs]", "E[μ_null]")

cenarios = [
 ("totals 150:900 Beta(8,6)", 150:900, Beta(8,6)),
 ("totals 150:900 Beta(2,2)", 150:900, Beta(2,2)),
 ("totals  20:120 Beta(8,6)",  20:120, Beta(8,6)),
 ("totals  50:400 Beta(8,6)",  50:400, Beta(8,6)),
 ("totals 300:400 Beta(8,6)", 300:400, Beta(8,6)),   # TSE: seções homogêneas
 ("totals 150:900 Beta(1.2,8)",150:900, Beta(1.2,8)),# vencedor minoritário
]

for (nome, tr, dist) in cenarios, m in (500, 2000)
    ps = Vector{Float64}(undef, R); Tobs = zeros(R); Tnull = zeros(R)
    Threads.@threads for r in 1:R
        rng = Xoshiro(hash((nome, m, r)))
        v, t = replica(rng, m, tr, dist)
        res = rozenas_test(v, t; B = B, rng = rng)
        ps[r] = res.pvalue; Tobs[r] = res.total_observed; Tnull[r] = mean(res.null_totals)
    end
    @printf("%-22s %-7d %-8.4f %-8.4f %-8.4f %-8.4f %-8.1f %-8.1f\n",
            nome, m, pvalue(ExactOneSampleKSTest(ps, Uniform(0,1))),
            mean(ps .< 0.01), mean(ps .< 0.05), mean(ps .< 0.10),
            mean(Tobs), mean(Tnull))
end
