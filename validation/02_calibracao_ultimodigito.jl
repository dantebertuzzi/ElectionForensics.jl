# Calibração sob H0 do teste de último/penúltimo dígito (Beber & Scacco 2012).
# H0 do teste: último dígito ~ Uniforme{0..9}. Isso só vale assintoticamente,
# quando o desvio-padrão da contagem é >> 10. Quantificamos onde quebra.
using ElectionForensics, Random, Distributions, Statistics, HypothesisTests, Printf

const R = 3000

println("="^78)
println("ÚLTIMO DÍGITO — erro tipo I sob eleição limpa (votos ~ Binomial(N, p))")
println("N = eleitorado da seção; sd(votos) = sqrt(N p (1-p)) governa a uniformidade")
println("="^78)
@printf("%-10s %-8s %-9s %-9s %-9s %-9s %-9s\n",
        "N seção", "n seções", "sd(votos)", "KS p", "rej@.01", "rej@.05", "rej@.10")

for N in (30, 60, 100, 200, 400, 1000, 5000), n in (1000,)
    ps = Vector{Float64}(undef, R)
    Threads.@threads for r in 1:R
        rng = Xoshiro(hash((:ld, N, n, r)))
        p = rand(rng, Beta(8, 6), n)
        x = [rand(rng, Binomial(N, p[i])) for i in 1:n]
        ps[r] = last_digit_test(x).pvalue
    end
    @printf("%-10d %-8d %-9.1f %-9.4f %-9.4f %-9.4f %-9.4f\n",
            N, n, sqrt(N*0.57*0.43), pvalue(ExactOneSampleKSTest(ps, Uniform(0,1))),
            mean(ps .< 0.01), mean(ps .< 0.05), mean(ps .< 0.10))
end

println()
println("Efeito do tamanho amostral n (N = 400, típico de seção do TSE):")
@printf("%-10s %-9s %-9s\n", "n seções", "KS p", "rej@.05")
for n in (500, 2000, 8000, 30000)
    ps = Vector{Float64}(undef, R)
    Threads.@threads for r in 1:R
        rng = Xoshiro(hash((:ldn, n, r)))
        p = rand(rng, Beta(8, 6), n)
        x = [rand(rng, Binomial(400, p[i])) for i in 1:n]
        ps[r] = last_digit_test(x).pvalue
    end
    @printf("%-10d %-9.4f %-9.4f\n", n, pvalue(ExactOneSampleKSTest(ps, Uniform(0,1))), mean(ps .< 0.05))
end

println()
println("min_value default = 10: seções pequenas passam pelo filtro.")
println("Mistura realista (N ~ 10..400, muitas seções pequenas):")
for n in (2000, 10000)
    ps = Vector{Float64}(undef, R)
    Threads.@threads for r in 1:R
        rng = Xoshiro(hash((:mix, n, r)))
        Ns = rand(rng, 10:400, n)
        p = rand(rng, Beta(8, 6), n)
        x = [rand(rng, Binomial(Ns[i], p[i])) for i in 1:n]
        ps[r] = last_digit_test(x).pvalue
    end
    @printf("  n=%-6d rej@.05 = %.4f   (nominal 0.05)\n", n, mean(ps .< 0.05))
end
