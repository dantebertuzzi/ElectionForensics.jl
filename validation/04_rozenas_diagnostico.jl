# Diagnósticos internos da nula do teste de Rozenas:
#  (a) a nula reproduz a variância verdadeira de T? (se inflada -> perda de poder)
#  (b) o clamp(share + h·z, 0, 1) cria átomos em 0 e 1 -> deflaciona a nula?
#  (c) sensibilidade a h e a B
using ElectionForensics, Random, Distributions, Statistics, Printf

const R = 600

function gen(rng, m, tr, dist)
    totals = rand(rng, tr, m)
    p = clamp.(rand(rng, dist, m), 1e-6, 1-1e-6)
    votes = [rand(rng, Binomial(totals[i], p[i])) for i in 1:m]
    votes, totals
end

println("="^88)
println("(a) sd verdadeiro de T (entre réplicas) vs sd estimado pela nula (dentro)")
println("="^88)
@printf("%-26s %-6s %-10s %-10s %-8s\n", "cenário", "m", "sd verdad.", "sd nula", "razão")
for (nome, tr, dist) in (("150:900 Beta(8,6)", 150:900, Beta(8,6)),
                         ("20:120 Beta(8,6)",   20:120, Beta(8,6)),
                         ("300:400 Beta(8,6)", 300:400, Beta(8,6)),
                         ("150:900 Beta(.5,.5)",150:900, Beta(0.5,0.5)))
    for m in (2000,)
        Tobs = zeros(R); sdn = zeros(R)
        Threads.@threads for r in 1:R
            rng = Xoshiro(hash((nome, m, r)))
            v, t = gen(rng, m, tr, dist)
            res = rozenas_test(v, t; B = 199, rng = rng)
            Tobs[r] = res.total_observed; sdn[r] = std(res.null_totals)
        end
        @printf("%-26s %-6d %-10.3f %-10.3f %-8.3f\n", nome, m, std(Tobs), mean(sdn), mean(sdn)/std(Tobs))
    end
end

println()
println("="^88)
println("(b) efeito do clamp: fração de p̃ colada em 0 ou 1 na geração da nula")
println("="^88)
for (nome, dist) in (("Beta(8,6)", Beta(8,6)), ("Beta(2,2)", Beta(2,2)),
                     ("Beta(.5,.5)", Beta(0.5,0.5)), ("Beta(1.2,8)", Beta(1.2,8)),
                     ("Beta(.3,.3) extremo", Beta(0.3,0.3)))
    rng = Xoshiro(42)
    v, t = gen(rng, 2000, 150:900, dist)
    shares = v ./ t
    h = ElectionForensics._silverman(shares)
    # replica exatamente o passo de jitter do pacote
    nclamp = 0; tot = 0
    for b in 1:200, i in eachindex(shares)
        p = shares[i] + h*randn(rng); tot += 1
        (p <= 0 || p >= 1) && (nclamp += 1)
    end
    @printf("  %-20s h=%.4f  p̃ clampado: %.3f%%\n", nome, h, 100*nclamp/tot)
end

println()
println("="^88)
println("(c) sensibilidade do p-valor a h e a B (mesmos dados; H0 verdadeiro)")
println("="^88)
rng = Xoshiro(7); v, t = gen(rng, 2000, 150:900, Beta(8,6))
hsil = ElectionForensics._silverman(v ./ t)
@printf("  h Silverman = %.4f\n", hsil)
for h in (hsil/8, hsil/4, hsil/2, hsil, 2hsil, 4hsil)
    ps = [rozenas_test(v, t; h = h, B = 999, rng = Xoshiro(1000+k)).pvalue for k in 1:5]
    mn = [mean(rozenas_test(v, t; h = h, B = 999, rng = Xoshiro(1000+k)).null_totals) for k in 1:3]
    @printf("  h=%.4f  p ∈ [%.3f, %.3f]  E[T|H0]≈%.1f  (T_obs=%d)\n",
            h, minimum(ps), maximum(ps), mean(mn),
            rozenas_test(v,t;h=h,B=99,rng=Xoshiro(1)).total_observed)
end
println()
println("  variação Monte Carlo do p-valor com B (mesmos dados, 10 sementes):")
for B in (99, 199, 999, 4999)
    ps = [rozenas_test(v, t; B = B, rng = Xoshiro(2000+k)).pvalue for k in 1:10]
    @printf("  B=%-6d p ∈ [%.4f, %.4f]  sd=%.4f\n", B, minimum(ps), maximum(ps), std(ps))
end
