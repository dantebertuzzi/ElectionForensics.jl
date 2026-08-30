# Procura um diagnóstico computável a partir dos dados que prediga quando o
# null uniforme (penúltimo dígito) e o null de Benford 2BL são aplicáveis.
# Candidato: S = log10(q95 / q05) = ordens de grandeza cobertas pelas contagens.
using ElectionForensics, Random, Distributions, Statistics, Printf

const R = 800
spread(x) = log10(quantile(Float64.(x), 0.95) / max(quantile(Float64.(x), 0.05), 1.0))

println("="^86)
println("S = log10(q95/q05) das contagens vs. erro tipo I (α = 0.05, m = 2000)")
println("="^86)
@printf("%-8s %-7s %-9s %-13s %-13s %-13s\n",
        "meanlog", "sdlog", "S médio", "rej penúlt.", "rej últ.", "rej 2BL")
for meanlog in (4.0, 6.0), sdlog in (0.1, 0.25, 0.5, 0.75, 1.0, 1.5, 2.0)
    pp = zeros(R); pl = zeros(R); pb = zeros(R); ss = zeros(R)
    Threads.@threads for r in 1:R
        rng = Xoshiro(hash((meanlog, sdlog, r)))
        totals = [max(30, round(Int, exp(meanlog + sdlog*randn(rng)))) for _ in 1:2000]
        p = rand(rng, Beta(8,6), 2000)
        v = [rand(rng, Binomial(totals[i], p[i])) for i in 1:2000]
        ss[r] = spread(v)
        pp[r] = any(≥(100), v) ? last_digit_test(v; position = :penultimate).pvalue : NaN
        pl[r] = last_digit_test(v; position = :last).pvalue
        pb[r] = benford_test(v; digit = 2).pvalue
    end
    @printf("%-8.1f %-7.2f %-9.2f %-13.4f %-13.4f %-13.4f\n",
            meanlog, sdlog, mean(ss), mean(filter(!isnan, pp) .< 0.05), mean(pl .< 0.05), mean(pb .< 0.05))
end
