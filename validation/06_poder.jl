# Curvas de poder: fraude injetada de intensidade crescente.
using ElectionForensics, Random, Distributions, Statistics, Printf

const R = 500

"Eleição limpa base."
function base(rng, m, tr)
    totals = rand(rng, tr, m)
    p = rand(rng, Beta(8, 6), m)
    votes = [rand(rng, Binomial(totals[i], p[i])) for i in 1:m]
    votes, totals
end

"Fraude 1: metas redondas — fração ε das seções forçada a k/d exato."
function fraude_metas!(rng, votes, totals, eps)
    alvos = (1//2, 3//5, 2//3, 7//10, 3//4, 4//5)
    for i in 1:round(Int, eps*length(votes))
        f = rand(rng, alvos)
        votes[i] = (totals[i] * numerator(f)) ÷ denominator(f)
    end
    votes
end

"Fraude 2: ballot stuffing — fração ε das seções recebe +δ% do total."
function fraude_stuffing!(rng, votes, totals, eps; delta = 0.15)
    for i in 1:round(Int, eps*length(votes))
        votes[i] = min(totals[i], votes[i] + round(Int, delta*totals[i]))
    end
    votes
end

"Fraude 3: arredondamento — contagens de ε das seções arredondadas à dezena."
function fraude_arred!(rng, votes, totals, eps)
    for i in 1:round(Int, eps*length(votes))
        votes[i] = min(totals[i], 10*round(Int, votes[i]/10))
    end
    votes
end

println("="^100)
println("PODER (α = 0.05, m = 2000, totals ∈ 150:900, B = 199, R = $R)")
println("="^100)
@printf("%-22s %-6s %-11s %-11s %-11s %-11s\n",
        "fraude", "ε", "rozenas", "últ.dígito", "penúlt.dig", "benford 2BL")

for (nome, inj) in (("metas redondas", fraude_metas!),
                    ("stuffing +15%",  fraude_stuffing!),
                    ("arredond. dezena", fraude_arred!))
    for eps in (0.0, 0.01, 0.02, 0.05, 0.10, 0.20)
        pr = zeros(R); pl = zeros(R); pp = zeros(R); pb = zeros(R)
        Threads.@threads for r in 1:R
            rng = Xoshiro(hash((nome, eps, r)))
            v, t = base(rng, 2000, 150:900)
            eps > 0 && inj(rng, v, t, eps)
            pr[r] = rozenas_test(v, t; B = 199, rng = rng).pvalue
            pl[r] = last_digit_test(v).pvalue
            pp[r] = last_digit_test(v; position = :penultimate).pvalue
            pb[r] = benford_test(v; digit = 2).pvalue
        end
        @printf("%-22s %-6.2f %-11.3f %-11.3f %-11.3f %-11.3f\n", nome, eps,
                mean(pr .< 0.05), mean(pl .< 0.05), mean(pp .< 0.05), mean(pb .< 0.05))
    end
    println()
end

println("Poder do teste de Rozenas vs. tamanho amostral (fraude 'metas redondas', ε = 0.02):")
@printf("  %-8s %-10s\n", "m", "poder")
for m in (250, 500, 1000, 2000, 4000)
    pr = zeros(R)
    Threads.@threads for r in 1:R
        rng = Xoshiro(hash((:m, m, r)))
        v, t = base(rng, m, 150:900)
        fraude_metas!(rng, v, t, 0.02)
        pr[r] = rozenas_test(v, t; B = 199, rng = rng).pvalue
    end
    @printf("  %-8d %-10.3f\n", m, mean(pr .< 0.05))
end
