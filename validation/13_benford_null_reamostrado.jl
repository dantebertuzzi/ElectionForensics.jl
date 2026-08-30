# Proposta: em vez de testar contra a lei de Benford (que contagens eleitorais
# limpas NÃO seguem), testar o 2º dígito contra um null REAMOSTRADO — kernel
# gaussiano em log10(contagem), banda >> escala do 2º dígito (~0.004 décadas)
# mas << amplitude dos dados. Pergunta respondida: "o 2º dígito é anômalo
# DADO o formato empírico da distribuição de contagens?"
using ElectionForensics, Random, Distributions, Statistics, Printf
const R = 500

second_dig(x) = (a = abs(x); while a >= 100; a ÷= 10; end; a % 10)

function chi2_2bl(v, expected)
    c = zeros(Int, 10); for x in v; c[second_dig(x)+1] += 1; end
    n = length(v); E = n .* expected
    sum((c .- E).^2 ./ E)
end

"null reamostrado: jitter gaussiano em log10, banda hlog décadas."
function benford2_resampled(v::Vector{Int}; B = 499, hlog = 0.05, rng = Random.default_rng())
    v = filter(>=(10), v); n = length(v)
    e = ElectionForensics.benford_expected(2)
    chi_obs = chi2_2bl(v, e)
    lv = log10.(Float64.(v))
    nulls = zeros(B)
    for b in 1:B
        vb = Vector{Int}(undef, 0); sizehint!(vb, n)
        while length(vb) < n
            y = round(Int, 10^(lv[rand(rng,1:n)] + hlog*randn(rng)))
            y >= 10 && push!(vb, y)
        end
        nulls[b] = chi2_2bl(vb, e)
    end
    (1 + count(>=(chi_obs), nulls)) / (B + 1)
end

dgps = [
 ("totals U(150,900) p~B(8,6)",   rng->(t=rand(rng,150:900,2000); p=rand(rng,Beta(8,6),2000); [rand(rng,Binomial(t[i],p[i])) for i in 1:2000])),
 ("totals U(100,9999) p~B(8,6)",  rng->(t=rand(rng,100:9999,2000); p=rand(rng,Beta(8,6),2000); [rand(rng,Binomial(t[i],p[i])) for i in 1:2000])),
 ("totals LN(6,1.2) p~B(8,6)",    rng->(t=[max(20,round(Int,exp(6+1.2randn(rng)))) for _ in 1:2000]; p=rand(rng,Beta(8,6),2000); [rand(rng,Binomial(t[i],p[i])) for i in 1:2000])),
 ("TSE-like: totals U(200,400)",  rng->(t=rand(rng,200:400,2000); p=rand(rng,Beta(8,6),2000); [rand(rng,Binomial(t[i],p[i])) for i in 1:2000])),
]
println("="^78)
println("2BL: erro tipo I (α=.05) — null de Benford (atual) vs null reamostrado")
println("="^78)
@printf("%-32s %-14s %-14s\n", "DGP", "2BL ATUAL", "2BL REAMOSTR.")
for (nome, gen) in dgps
    a=zeros(R); b=zeros(R)
    Threads.@threads for r in 1:R
        rng = Xoshiro(hash((nome,r))); v = gen(rng)
        a[r]=benford_test(v; digit=2).pvalue
        b[r]=benford2_resampled(v; B=299, rng=rng)
    end
    @printf("%-32s %-14.4f %-14.4f\n", nome, mean(a.<0.05), mean(b.<0.05))
end
println()
println("PODER: contagens de ε das seções arredondadas à dezena (TSE-like):")
@printf("%-8s %-14s %-14s\n", "ε", "2BL ATUAL", "2BL REAMOSTR.")
for eps in (0.0, 0.05, 0.15, 0.30)
    a=zeros(R); b=zeros(R)
    Threads.@threads for r in 1:R
        rng = Xoshiro(hash((:p,eps,r)))
        t=rand(rng,200:400,2000); p=rand(rng,Beta(8,6),2000)
        v=[rand(rng,Binomial(t[i],p[i])) for i in 1:2000]
        for i in 1:round(Int,eps*2000); v[i]=10*round(Int,v[i]/10); end
        a[r]=benford_test(v; digit=2).pvalue
        b[r]=benford2_resampled(v; B=299, rng=rng)
    end
    @printf("%-8.2f %-14.3f %-14.3f\n", eps, mean(a.<0.05), mean(b.<0.05))
end

println()
println("PODER contra manipulação que ATINGE o 2º dígito (arredondar à centena):")
@printf("%-8s %-14s %-14s\n", "ε", "2BL ATUAL", "2BL REAMOSTR.")
for eps in (0.0, 0.02, 0.05, 0.10, 0.20)
    a=zeros(R); b=zeros(R)
    Threads.@threads for r in 1:R
        rng = Xoshiro(hash((:c,eps,r)))
        t=rand(rng,200:900,2000); p=rand(rng,Beta(8,6),2000)
        v=[rand(rng,Binomial(t[i],p[i])) for i in 1:2000]
        for i in 1:round(Int,eps*2000); v[i]=max(10, 100*round(Int,v[i]/100)); end
        a[r]=benford_test(v; digit=2).pvalue
        b[r]=benford2_resampled(v; B=299, rng=rng)
    end
    @printf("%-8.2f %-14.3f %-14.3f\n", eps, mean(a.<0.05), mean(b.<0.05))
end
