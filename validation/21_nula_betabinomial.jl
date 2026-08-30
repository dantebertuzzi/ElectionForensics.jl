# A nula por mistura Beta-Binomial (modelo do `spikes`, de Rozenas) contra a
# nula por kernel: calibração sob H0 e poder sob fraude injetada.
using ElectionForensics, Random, Distributions, Statistics, Printf
const R = 400
const NULAS = (:kernel, :betabinomial)

println("="^96)
println("ERRO TIPO I sob H0 (α = 0,05, m = 1500, B = 199).  Alvo: 0,050")
println("="^96)
@printf("%-26s %-14s %-14s\n", "distribuição de shares", "kernel", "beta-binomial")
for (nome, dist) in (("Beta(8,6)",    Beta(8,6)),
                     ("Beta(2,2)",    Beta(2,2)),
                     ("Beta(1,1)",    Beta(1,1)),
                     ("Beta(.5,.5)",  Beta(0.5,0.5)),
                     ("Beta(.3,.3)",  Beta(0.3,0.3)),
                     ("Beta(1.2,8)",  Beta(1.2,8)),
                     ("bimodal 2 comp", MixtureModel([Beta(20,30), Beta(45,15)], [0.5,0.5])))
    res = Dict(k => zeros(R) for k in NULAS)
    Threads.@threads for r in 1:R
        rng = Xoshiro(hash((nome, r)))
        t = rand(rng, 150:900, 1500)
        p = clamp.(rand(rng, dist, 1500), 1e-9, 1-1e-9)
        v = [rand(rng, Binomial(t[i], p[i])) for i in 1:1500]
        for k in NULAS
            res[k][r] = rozenas_test(v, t; B = 199, null = k,
                                     rng = Xoshiro(hash((nome, r, k)))).pvalue
        end
    end
    @printf("%-26s %-14.4f %-14.4f\n", nome,
            mean(res[:kernel] .< 0.05), mean(res[:betabinomial] .< 0.05))
end

println()
println("="^96)
println("PODER: metas redondas em ε das seções (m = 1500)")
println("="^96)
@printf("%-16s %-8s %-14s %-14s\n", "shares", "ε", "kernel", "beta-binomial")
for (nome, dist) in (("Beta(8,6)", Beta(8,6)), ("Beta(.5,.5)", Beta(0.5,0.5)))
    for eps in (0.01, 0.02, 0.05)
        res = Dict(k => zeros(R) for k in NULAS)
        Threads.@threads for r in 1:R
            rng = Xoshiro(hash((:pw, nome, eps, r)))
            t = rand(rng, 150:900, 1500)
            p = clamp.(rand(rng, dist, 1500), 1e-9, 1-1e-9)
            v = [rand(rng, Binomial(t[i], p[i])) for i in 1:1500]
            alvos = (1//2, 3//5, 2//3, 7//10, 3//4, 4//5)
            for i in 1:round(Int, eps*1500)
                f = rand(rng, alvos); v[i] = (t[i]*numerator(f)) ÷ denominator(f)
            end
            for k in NULAS
                res[k][r] = rozenas_test(v, t; B = 199, null = k,
                                         rng = Xoshiro(hash((:pw,nome,eps,r,k)))).pvalue
            end
        end
        @printf("%-16s %-8.2f %-14.3f %-14.3f\n", nome, eps,
                mean(res[:kernel] .< 0.05), mean(res[:betabinomial] .< 0.05))
    end
end
