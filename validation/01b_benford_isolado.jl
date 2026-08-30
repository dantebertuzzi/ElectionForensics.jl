# Isola: (i) a maquinaria χ² está correta? (ii) a inflação vem do DGP?
using ElectionForensics, Random, Distributions, Statistics, HypothesisTests, Printf

const R = 5000

# ── (i) construir inteiros cujo 2º dígito vem EXATAMENTE da lei de Benford
#       (amostrando o dígito e completando com ruído uniforme)
e2 = ElectionForensics.benford_expected(2)
e1 = ElectionForensics.benford_expected(1)
println("── (i) dígitos amostrados diretamente da lei de Benford ──")
for (d, ev) in ((1, e1), (2, e2)), n in (500, 2000)
    ps = Vector{Float64}(undef, R)
    Threads.@threads for r in 1:R
        rng = Xoshiro(hash((:pure, d, n, r)))
        cat = Categorical(ev)
        x = if d == 1
            [rand(rng, cat) * 100 + rand(rng, 0:99) for _ in 1:n]       # d1 = dígito
        else
            [rand(rng, 1:9) * 100 + (rand(rng, cat) - 1) * 10 + rand(rng, 0:9) for _ in 1:n]
        end
        ps[r] = benford_test(x; digit = d).pvalue
    end
    @printf("  d=%d n=%5d | KS p=%.4f | α=.05 → %.4f (ep=%.4f)\n",
            d, n, pvalue(ExactOneSampleKSTest(ps, Uniform(0,1))),
            mean(ps .< 0.05), sqrt(0.05*0.95/R))
end

# ── (ii) sensibilidade da 2BL à dispersão dos tamanhos de seção ──────
println()
println("── (ii) DGP eleitoral limpo: taxa de erro tipo I da 2BL (α = 0.05) ──")
println("   sdlog controla quantas ordens de grandeza os tamanhos de seção cobrem")
@printf("  %-8s %-8s %-10s %-10s\n", "sdlog", "n", "rej@.05", "MAD médio")
for sdlog in (0.2, 0.5, 1.0, 1.5, 2.0), n in (2000,)
    ps = Vector{Float64}(undef, R ÷ 5); mads = similar(ps)
    Threads.@threads for r in 1:(R ÷ 5)
        rng = Xoshiro(hash((:sd, sdlog, n, r)))
        totals = [max(20, round(Int, exp(6.0 + sdlog*randn(rng)))) for _ in 1:n]
        p = rand(rng, Beta(8, 6), n)
        x = [rand(rng, Binomial(totals[i], p[i])) for i in 1:n]
        res = benford_test(x; digit = 2)
        ps[r] = res.pvalue; mads[r] = res.mad
    end
    @printf("  %-8.1f %-8d %-10.4f %-10.5f\n", sdlog, n, mean(ps .< 0.05), mean(mads))
end
