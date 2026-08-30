# `calibration_check` acerta? Comparamos o que ele PREVÊ com a taxa de erro
# tipo I MEDIDA por simulação independente, nos mesmos regimes.
using ElectionForensics, Random, Distributions, Statistics, Printf
const R = 600
println("="^104)
println("VALIDAÇÃO DO AUTO-DIAGNÓSTICO: previsto vs. medido (α = 0,05)")
println("="^104)
@printf("%-26s %-11s %-11s %-11s %-11s %-14s\n","regime","últ MEDIDO","últ PREVISTO",
        "pen MEDIDO","pen PREVISTO","veredito últ")
for (tr,dist,rot) in [(150:900, Beta(8,6),     "B(8,6) tot 150:900"),
                      (150:900, Beta(0.3,0.3), "B(.3,.3) tot 150:900"),
                      (20:400,  Beta(0.3,0.3), "B(.3,.3) tot 20:400"),
                      (10:120,  Beta(8,6),     "B(8,6) tot 10:120"),
                      (1000:9000, Beta(8,6),   "B(8,6) tot 1000:9000")]
    # medido: R eleições limpas independentes
    pl=zeros(R); pp=zeros(R)
    Threads.@threads for r in 1:R
        rng=Xoshiro(hash((rot,r)))
        t=rand(rng,tr,2500); p=clamp.(rand(rng,dist,2500),1e-9,1-1e-9)
        v=[rand(rng,Binomial(t[i],p[i])) for i in 1:2500]
        pl[r]=last_digit_test(v; warn=false).pvalue
        pp[r]= any(≥(100),v) ? last_digit_test(v; position=:penultimate, warn=false).pvalue : NaN
    end
    medido_l=mean(pl.<0.05); medido_p=mean(filter(!isnan,pp).<0.05)
    # previsto: uma única amostra observada, diagnosticada pelo pacote
    rng=Xoshiro(hash((rot,:obs)))
    t=rand(rng,tr,2500); p=clamp.(rand(rng,dist,2500),1e-9,1-1e-9)
    v=[rand(rng,Binomial(t[i],p[i])) for i in 1:2500]
    cal = calibration_check(v,t; R=250, B=99, tests=[:last_digit,:penultimate], rng=rng)
    @printf("%-26s %-11.3f %-11.3f %-11.3f %-11.3f %-14s\n", rot,
            medido_l, cal.rejection_rate[1], medido_p, cal.rejection_rate[2],
            String(cal.verdict[1]))
end
