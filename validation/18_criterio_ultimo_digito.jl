# Existe um diagnóstico calculável a partir dos dados que prediga quando o null
# uniforme do ÚLTIMO dígito falha? Varremos regimes e correlacionamos o erro
# tipo I com candidatos a estatística de aplicabilidade.
using ElectionForensics, Random, Distributions, Statistics, Printf
const R = 800
mv = 10

# candidatos a diagnóstico, calculados sobre as contagens RETIDAS
frac_pequenas(v) = count(<(50), v) / length(v)       # proporção abaixo de 50
q05(v)           = quantile(Float64.(v), 0.05)       # 5º percentil
sd_rel(v)        = std(v) / 10                       # dispersão em unidades de 10

println("="^112)
println("Erro tipo I do último dígito (mv=10) vs. candidatos a diagnóstico, α = 0,05, m = 3000")
println("="^112)
@printf("%-30s %-9s %-11s %-9s %-9s\n","regime","rej@.05","%v<50","q05(v)","sd/10")
linhas = []
for (tr, dist, rot) in [(150:900, Beta(8,6),    "B(8,6)  tot 150:900"),
                        (150:900, Beta(2,2),    "B(2,2)  tot 150:900"),
                        (150:900, Beta(1,1),    "B(1,1)  tot 150:900"),
                        (150:900, Beta(0.5,0.5),"B(.5,.5) tot 150:900"),
                        (150:900, Beta(0.3,0.3),"B(.3,.3) tot 150:900"),
                        (20:400,  Beta(8,6),    "B(8,6)  tot 20:400"),
                        (20:400,  Beta(0.3,0.3),"B(.3,.3) tot 20:400"),
                        (50:200,  Beta(8,6),    "B(8,6)  tot 50:200"),
                        (300:400, Beta(8,6),    "B(8,6)  tot 300:400"),
                        (10:120,  Beta(8,6),    "B(8,6)  tot 10:120"),
                        (1000:9000, Beta(8,6),  "B(8,6)  tot 1000:9000")]
    ps=zeros(R); f=zeros(R); q=zeros(R); s=zeros(R)
    Threads.@threads for r in 1:R
        rng=Xoshiro(hash((rot,r)))
        t=rand(rng,tr,3000); p=clamp.(rand(rng,dist,3000),1e-9,1-1e-9)
        v=[rand(rng,Binomial(t[i],p[i])) for i in 1:3000]
        vv=filter(≥(mv),v)
        if length(vv) < 200; ps[r]=NaN; continue; end
        ps[r]=last_digit_test(v; min_value=mv, warn=false).pvalue
        f[r]=frac_pequenas(vv); q[r]=q05(vv); s[r]=sd_rel(vv)
    end
    ok=.!isnan.(ps)
    push!(linhas,(rot, mean(ps[ok].<0.05), mean(f[ok]), mean(q[ok]), mean(s[ok])))
    @printf("%-30s %-9.4f %-11.3f %-9.0f %-9.1f\n", rot, linhas[end][2], linhas[end][3], linhas[end][4], linhas[end][5])
end

println()
println("Correlação de Spearman entre o erro tipo I e cada candidato:")
sp(x,y)=cor(sortperm(sortperm(x)), sortperm(sortperm(y)))
rej=[l[2] for l in linhas]
@printf("  %%v<50  : %+.3f\n", sp(rej,[l[3] for l in linhas]))
@printf("  q05(v) : %+.3f\n", sp(rej,[l[4] for l in linhas]))
@printf("  sd/10  : %+.3f\n", sp(rej,[l[5] for l in linhas]))
