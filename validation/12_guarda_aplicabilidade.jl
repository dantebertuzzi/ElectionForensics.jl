# Um guarda computável a partir dos dados salva a 2BL?
# S = log10(q95/q05) das contagens retidas. Testamos o critério S ≥ 1.2
# em DGPs variados (totals uniformes, lognormais, Beta diferentes).
using ElectionForensics, Random, Distributions, Statistics, Printf
const R = 600
S(v) = log10(quantile(Float64.(v), 0.95) / max(quantile(Float64.(v), 0.05), 1.0))

println("="^80)
println("2BL: erro tipo I (α=.05) vs S = log10(q95/q05), vários DGPs, m = 2000")
println("="^80)
@printf("%-34s %-8s %-10s %-10s\n", "DGP", "S médio", "rej 2BL", "rej 1BL")
dgps = [
 ("totals U(150,900) p~B(8,6)",   rng->(t=rand(rng,150:900,2000); p=rand(rng,Beta(8,6),2000); [rand(rng,Binomial(t[i],p[i])) for i in 1:2000])),
 ("totals U(100,9999) p~B(8,6)",  rng->(t=rand(rng,100:9999,2000); p=rand(rng,Beta(8,6),2000); [rand(rng,Binomial(t[i],p[i])) for i in 1:2000])),
 ("totals U(10,99999) p~B(8,6)",  rng->(t=rand(rng,10:99999,2000); p=rand(rng,Beta(8,6),2000); [rand(rng,Binomial(t[i],p[i])) for i in 1:2000])),
 ("totals LN(6,0.5) p~B(8,6)",    rng->(t=[max(20,round(Int,exp(6+0.5randn(rng)))) for _ in 1:2000]; p=rand(rng,Beta(8,6),2000); [rand(rng,Binomial(t[i],p[i])) for i in 1:2000])),
 ("totals LN(6,1.2) p~B(8,6)",    rng->(t=[max(20,round(Int,exp(6+1.2randn(rng)))) for _ in 1:2000]; p=rand(rng,Beta(8,6),2000); [rand(rng,Binomial(t[i],p[i])) for i in 1:2000])),
 ("totals LN(7,1.8) p~B(1.2,8)",  rng->(t=[max(20,round(Int,exp(7+1.8randn(rng)))) for _ in 1:2000]; p=rand(rng,Beta(1.2,8),2000); [rand(rng,Binomial(t[i],p[i])) for i in 1:2000])),
 ("Pareto(1) ×1000 (cauda pesada)",rng->(t=[max(20,round(Int,1000*rand(rng,Pareto(1.0)))) for _ in 1:2000]; p=rand(rng,Beta(8,6),2000); [rand(rng,Binomial(t[i],p[i])) for i in 1:2000])),
 ("TSE-like: totals U(200,400)",  rng->(t=rand(rng,200:400,2000); p=rand(rng,Beta(8,6),2000); [rand(rng,Binomial(t[i],p[i])) for i in 1:2000])),
]
for (nome, gen) in dgps
    p2=zeros(R); p1=zeros(R); ss=zeros(R)
    Threads.@threads for r in 1:R
        rng = Xoshiro(hash((nome,r))); v = gen(rng)
        ss[r]=S(filter(>=(10), v))
        p2[r]=benford_test(v; digit=2).pvalue
        p1[r]=benford_test(v; digit=1).pvalue
    end
    @printf("%-34s %-8.2f %-10.4f %-10.4f\n", nome, mean(ss), mean(p2.<0.05), mean(p1.<0.05))
end
println()
println("⇒ guarda proposto: avisar quando S < 1.2 (2BL) — abaixo disso o null")
println("  de Benford não descreve contagens eleitorais limpas.")
