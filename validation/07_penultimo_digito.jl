# O penúltimo dígito só é ~uniforme se a distribuição das contagens for suave
# numa escala de 100. Com sd(votos) << 100 isso falha e o χ² rejeita H0 em
# dados limpos. Quantificamos o erro tipo I em função de sd(votos).
using ElectionForensics, Random, Distributions, Statistics, Printf

const R = 2000

println("="^92)
println("PENÚLTIMO DÍGITO — erro tipo I sob eleição limpa (min_value default = 100)")
println("="^92)
@printf("%-24s %-10s %-10s %-9s %-9s %-9s\n",
        "cenário", "sd(votos)", "amplitude", "rej@.01", "rej@.05", "rej@.10")

cenarios = [
  ("totals 150:900 B(8,6)", 150:900, Beta(8,6)),
  ("totals 300:400 B(8,6)", 300:400, Beta(8,6)),
  ("totals 100:200 B(8,6)", 100:200, Beta(8,6)),
  ("totals 100:9999 B(8,6)", 100:9999, Beta(8,6)),
  ("totals 1000:99999 B(8,6)", 1000:99999, Beta(8,6)),
  ("totals 150:900 B(1,1)", 150:900, Beta(1,1)),
]
for (nome, tr, dist) in cenarios
    ps = Vector{Float64}(undef, R); sds = zeros(R); amp = zeros(R)
    Threads.@threads for r in 1:R
        rng = Xoshiro(hash((nome, r)))
        totals = rand(rng, tr, 2000)
        p = rand(rng, dist, 2000)
        v = [rand(rng, Binomial(totals[i], p[i])) for i in 1:2000]
        ps[r] = last_digit_test(v; position = :penultimate).pvalue
        sds[r] = std(v); amp[r] = maximum(v) - minimum(v)
    end
    @printf("%-24s %-10.1f %-10.0f %-9.4f %-9.4f %-9.4f\n",
            nome, mean(sds), mean(amp), mean(ps .< 0.01), mean(ps .< 0.05), mean(ps .< 0.10))
end

println()
println("Mesmo diagnóstico para o ÚLTIMO dígito (controle — deve ficar em α):")
for (nome, tr, dist) in cenarios
    ps = Vector{Float64}(undef, R)
    Threads.@threads for r in 1:R
        rng = Xoshiro(hash((nome, :l, r)))
        totals = rand(rng, tr, 2000); p = rand(rng, dist, 2000)
        v = [rand(rng, Binomial(totals[i], p[i])) for i in 1:2000]
        ps[r] = last_digit_test(v; position = :last).pvalue
    end
    @printf("%-24s rej@.05 = %.4f\n", nome, mean(ps .< 0.05))
end

println()
println("Distribuição média do penúltimo dígito, cenário 150:900 B(8,6) (esp. 10%% cada):")
acc = zeros(10)
for r in 1:400
    rng = Xoshiro(hash((:dist, r)))
    totals = rand(rng, 150:900, 2000); p = rand(rng, Beta(8,6), 2000)
    v = [rand(rng, Binomial(totals[i], p[i])) for i in 1:2000]
    acc .+= last_digit_test(v; position = :penultimate).observed
end
acc ./= 400
for d in 0:9
    @printf("  %d: %.2f%%  %s\n", d, 100acc[d+1], "#"^round(Int, 400acc[d+1]))
end
