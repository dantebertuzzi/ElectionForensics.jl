# O `show` de RozenasResult imprime "frações com maior excesso" ordenadas por z.
# Esses z são MÁXIMOS de ~31 estatísticas correlacionadas: sob H0 já são grandes.
# Sem correção de multiplicidade, o usuário lê ruído como evidência.
using ElectionForensics, Random, Distributions, Statistics, Printf

const R = 1000
fr = coarse_fractions(10)
println("nº de frações coarse com max_denom = 10: ", length(fr))
println()
maxz = zeros(R); n_acima2 = zeros(Int, R)
Threads.@threads for r in 1:R
    rng = Xoshiro(hash((:z, r)))
    totals = rand(rng, 150:900, 2000)
    p = rand(rng, Beta(8,6), 2000)
    v = [rand(rng, Binomial(totals[i], p[i])) for i in 1:2000]
    res = rozenas_test(v, totals; B = 499, rng = rng)
    zs = [res.null_sd[j] > 0 ? (res.observed[j]-res.null_mean[j])/res.null_sd[j] : 0.0
          for j in eachindex(res.fractions)]
    maxz[r] = maximum(zs); n_acima2[r] = count(>(2.0), zs)
end
println("SOB H0 (eleições limpas), distribuição do maior z exibido:")
@printf("  mediana = %.2f   q90 = %.2f   q99 = %.2f   máx = %.2f\n",
        median(maxz), quantile(maxz,0.9), quantile(maxz,0.99), maximum(maxz))
@printf("  P(algum z > 2) = %.3f    P(algum z > 3) = %.3f\n",
        mean(maxz .> 2), mean(maxz .> 3))
@printf("  nº médio de frações com z > 2 por eleição limpa: %.2f\n", mean(n_acima2))
println()
println("⇒ ver 'z = 2,5' numa fração isolada NÃO é evidência: acontece em ",
        round(100*mean(maxz .> 2.5)), "% das eleições limpas.")

# p-valor de Rozenas em paralelo (controle: o teste agregado está calibrado)
