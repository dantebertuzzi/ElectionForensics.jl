# Exporta datasets fixos + as estatísticas do pacote, para conferência externa.
using ElectionForensics, Random, Distributions, Printf
mkpath(joinpath(@__DIR__, "xref"))
rng = Xoshiro(20260829)
datasets = Dict{String,Vector{Int}}(
    "loguniforme" => [round(Int, 10.0^(1 + 5rand(rng))) for _ in 1:5000],
    "binomial"    => [rand(rng, Binomial(rand(rng,150:900), rand(rng, Beta(8,6)))) for _ in 1:5000],
    "arredondado" => [10*rand(rng, 10:999) for _ in 1:5000],
)
open(joinpath(@__DIR__, "xref", "dados.csv"), "w") do io
    println(io, "dataset,valor")
    for (k, v) in sort(collect(datasets), by=first), x in v
        println(io, k, ",", x)
    end
end
open(joinpath(@__DIR__, "xref", "julia_resultados.csv"), "w") do io
    println(io, "dataset,teste,chi2,df,pvalue,mad")
    for (k, v) in sort(collect(datasets), by=first)
        for d in (1, 2)
            r = benford_test(v; digit=d)
            @printf(io, "%s,benford%d,%.12g,%d,%.12g,%.12g\n", k, d, r.chi2, r.df, r.pvalue, r.mad)
        end
        for pos in (:last, :penultimate)
            r = last_digit_test(v; position=pos)
            @printf(io, "%s,%s,%.12g,%d,%.12g,NaN\n", k, pos, r.chi2, r.df, r.pvalue)
        end
    end
end
println("exportado")
