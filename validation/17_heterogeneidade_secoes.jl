# Hipótese: o teste de Rozenas é anti-conservador quando os tamanhos de seção
# são MUITO heterogêneos (de N=2 a N=5000), porque a nula jitteia em torno de
# shares observados que, nas seções minúsculas, são puro ruído.
# Testamos com a distribuição EMPÍRICA de totais dos dados americanos reais.
using ElectionForensics, Random, Distributions, Statistics, Printf
function ler(path; cv=1, ct=2)
    a=Int[]; b=Int[]
    for (k,l) in enumerate(eachline(path))
        k==1 && continue; c=split(l,',')
        x=tryparse(Int,c[cv]); y=tryparse(Int,c[ct])
        (x===nothing||y===nothing) && continue
        (y>0 && 0≤x≤y) || continue; push!(a,x); push!(b,y)
    end; a,b
end
D=ARGS[1]; R=800
println("="^100)
println("ERRO TIPO I sob H0 com totais REAIS (bootstrap dos totals de cada UF)")
println("shares latentes ~ Beta suave ⇒ H0 verdadeiro por construção. α = 0,05")
println("="^100)
@printf("%-6s %-8s %-9s %-9s %-11s %-11s\n","UF","n","min(N)","%N≤100","rej TODOS","rej N≥100")
for st in ("ak","sd","me","mt","nh","vt")
    _, tot = ler(joinpath(D,"us2020_$st.csv"))
    m = length(tot)
    p_all = zeros(R); p_big = zeros(R)
    grandes = filter(≥(100), tot)
    Threads.@threads for r in 1:R
        rng = Xoshiro(hash((st,r)))
        # reamostra os totais reais; shares latentes de uma Beta suave
        T  = [tot[rand(rng,1:m)] for _ in 1:m]
        pp = rand(rng, Beta(8,6), m)
        V  = [rand(rng, Binomial(T[i], pp[i])) for i in 1:m]
        p_all[r] = rozenas_test(V, T; B=199, rng=rng).pvalue
        mg = length(grandes)
        if mg > 50
            T2 = [grandes[rand(rng,1:mg)] for _ in 1:mg]
            p2 = rand(rng, Beta(8,6), mg)
            V2 = [rand(rng, Binomial(T2[i], p2[i])) for i in 1:mg]
            p_big[r] = rozenas_test(V2, T2; B=199, rng=rng).pvalue
        else
            p_big[r] = NaN
        end
    end
    @printf("%-6s %-8d %-9d %-9.0f %-11.4f %-11s\n", uppercase(st), m, minimum(tot),
            100count(≤(100),tot)/m, mean(p_all.<0.05),
            all(isnan,p_big) ? "—" : @sprintf("%.4f", mean(filter(!isnan,p_big).<0.05)))
end
