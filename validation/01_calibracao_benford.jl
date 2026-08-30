# Calibração sob H0 do teste de Benford.
#   (a) DGP log-uniforme -> satisfaz Benford exatamente (valida a maquinaria χ²)
#   (b) DGP eleitoral limpo -> votos ~ Binomial(total, p); é o caso de uso real
using ElectionForensics, Random, Distributions, Statistics, HypothesisTests, Printf

"p-valores devem ser U(0,1) sob H0; testamos com KS."
function ks_unif(ps)
    t = ExactOneSampleKSTest(ps, Uniform(0, 1))
    return pvalue(t)
end

function taxas(ps)
    (a01 = mean(ps .< 0.01), a05 = mean(ps .< 0.05), a10 = mean(ps .< 0.10))
end

# ── (a) DGP log-uniforme: Benford exato ──────────────────────────────
function dgp_loguniforme(rng, n)
    [round(Int, 10.0^(1 + 5rand(rng))) for _ in 1:n]
end

# ── (b) DGP eleitoral limpo: seções lognormais, p ~ Beta ─────────────
function dgp_eleicao(rng, n; meanlog = 6.0, sdlog = 0.5, a = 8.0, b = 6.0)
    totals = [max(20, round(Int, exp(meanlog + sdlog * randn(rng)))) for _ in 1:n]
    p = rand(rng, Beta(a, b), n)
    [rand(rng, Binomial(totals[i], p[i])) for i in 1:n]
end

const R = 2000

println("="^70)
println("BENFORD — calibração sob H0 (R = $R réplicas)")
println("="^70)

for (nome, dgp) in (("log-uniforme (Benford exato)", dgp_loguniforme),
                    ("eleição limpa (Binomial)", dgp_eleicao))
    for digit in (1, 2), n in (500, 2000)
        ps = Vector{Float64}(undef, R)
        Threads.@threads for r in 1:R
            rng = Xoshiro(hash((nome, digit, n, r)))   # semente por réplica: reprodutível
            x = dgp(rng, n)
            ps[r] = benford_test(x; digit = digit).pvalue
        end
        t = taxas(ps)
        @printf("%-30s d=%d n=%5d | KS p=%.4f | α=.01→%.3f  α=.05→%.3f  α=.10→%.3f\n",
                nome, digit, n, ks_unif(ps), t.a01, t.a05, t.a10)
    end
end
