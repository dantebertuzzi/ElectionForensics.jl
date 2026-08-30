# Casos-limite: n pequeno, zeros, empates, eleitorado uniforme, dados
# já arredondados, missing/NaN, tipos inteiros exóticos, overflow.
using ElectionForensics, Random, Distributions, Statistics, Printf

res(lbl, f) = println(rpad(lbl, 52), " → ", try
        v = f(); v isa Nothing ? "OK" : string(v)
    catch e
        "⚠ " * string(nameof(typeof(e))) * ": " *
        first(split(sprint(showerror, e), "\n"))[1:min(end,70)]
    end)

println("="^100); println("CASOS-LIMITE"); println("="^100)

res("benford_test n=1 (digit=2)",            () -> benford_test([10]; digit=2).pvalue)
res("benford_test todos iguais (x=[50]*500)",() -> round(benford_test(fill(50,500);digit=2).pvalue, digits=6))
res("benford_test com zeros (metade)",       () -> (r=benford_test([zeros(Int,500); rand(10:99999,500)];digit=1); (n=r.n, exc=r.n_excluded)))
res("benford_test valores negativos",        () -> (r=benford_test([-5,-100,12,345,6789];digit=1); (n=r.n, exc=r.n_excluded)))
res("last_digit_test n=1",                   () -> last_digit_test([10]).pvalue)
res("last_digit_test todos = 100",           () -> last_digit_test(fill(100,500)).pvalue)
res("last_digit_test min_value negativo",    () -> (r=last_digit_test([-13,-27,40,55]; min_value=-100); r.counts))
res("rozenas m=1",                           () -> rozenas_test([50],[100]; B=99).pvalue)
res("rozenas todos empatados 1//2",          () -> (r=rozenas_test(fill(50,500), fill(100,500); B=99); (T=r.total_observed, p=r.pvalue, h=r.h)))
res("rozenas todos share=0",                 () -> (r=rozenas_test(zeros(Int,500), fill(100,500); B=99); (T=r.total_observed, p=r.pvalue, h=r.h)))
res("rozenas todos share=1 (unânime)",       () -> (r=rozenas_test(fill(100,500), fill(100,500); B=99); (T=r.total_observed, p=r.pvalue, h=r.h)))
res("rozenas totals=1 (eleitorado unitário)",() -> (r=rozenas_test(rand(0:1,500), ones(Int,500); B=99); (T=r.total_observed, p=r.pvalue)))
res("rozenas max_denom grande (50)",         () -> length(coarse_fractions(50)))
res("rozenas totals gigantes (10^9)",        () -> (r=rozenas_test([5*10^8], [10^9]; B=99); r.total_observed))

println()
println("── missing / NaN / tipos ──")
res("benford_test Vector{Union{Missing,Int}}", () -> benford_test(Union{Missing,Int}[10,20,missing]))
res("benford_test Vector{Float64}",            () -> benford_test([10.0,20.0,30.0]))
res("benford_test Vector{Int8}",               () -> benford_test(Int8[10,20,30,40,50,60];digit=2).n)
res("benford_test Vector{BigInt}",             () -> benford_test(BigInt[10,20,30];digit=2).n)
res("rozenas Int16 com totals grandes",        () -> rozenas_test(Int16[100], Int16[300]; B=99).total_observed)
res("last_digit_test Bool",                    () -> last_digit_test([true,false]))

println()
println("── dados já arredondados (comuns em boletins agregados) ──")
rng = Xoshiro(3)
v = [10*rand(rng, 10:99) for _ in 1:2000]
r = last_digit_test(v)
@printf("  todas contagens múltiplas de 10: último dígito p = %.3e, freq{0,5} = %.3f\n", r.pvalue, r.freq_0_5)
r2 = benford_test(v; digit=2)
@printf("  mesmas contagens, 2BL: p = %.4f, conformidade = %s\n", r2.pvalue, r2.conformity)

println()
println("── determinismo / reprodutibilidade ──")
tot = rand(Xoshiro(9), 150:900, 800); vv = [rand(Xoshiro(hash(i)), Binomial(tot[i],0.55)) for i in eachindex(tot)]
p1 = rozenas_test(vv, tot; B=299, rng=Xoshiro(123)).pvalue
p2 = rozenas_test(vv, tot; B=299, rng=Xoshiro(123)).pvalue
println("  mesmo rng ⇒ mesmo p: ", p1 == p2, "  (p = ", p1, ")")
println("  nthreads = ", Threads.nthreads(), " (pacote é single-thread; sem dependência de threads)")
