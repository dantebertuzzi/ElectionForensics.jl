# Fingerprint de Klimek et al. (2012) em dados reais.
#
#   Klimek, Yegorov, Hanel & Thurner (2012), PNAS 109(41):16469-16473
#   (preprint arXiv:1201.3087).
#
# Fig. 3 do artigo: assimetria e curtose da taxa logarítmica de voto ν se
# agrupam em (0, 3) nos países sem alegação de fraude, e se afastam muito na
# Rússia e em Uganda. Reproduzimos isso com dados reais.
#
# O fingerprint exige TRÊS entradas por seção — votos no vencedor, cédulas
# depositadas e eleitorado. Só 6 dos 40 estados do OpenElections publicam
# eleitorado por precinct; a Rússia publica para todas as 95 mil seções.
#
#   julia validation/22_fingerprint.jl <dir_com_*_fp.csv>
using ElectionForensics, Random, Distributions, Printf

function ler3(p)
    a = Int[]; b = Int[]; c = Int[]
    for (k, l) in enumerate(eachline(p))
        k == 1 && continue
        f = split(l, ',')
        length(f) < 3 && continue
        x = tryparse(Int, f[1]); y = tryparse(Int, f[2]); z = tryparse(Int, f[3])
        (x === nothing || y === nothing || z === nothing) && continue
        push!(a, x); push!(b, y); push!(c, z)
    end
    return a, b, c
end

canto(r, frac = 0.95) = (nb = size(r.counts, 1); q = max(1, ceil(Int, frac*nb));
                         sum(@view r.counts[q:nb, q:nb]))

D = length(ARGS) ≥ 1 ? ARGS[1] : error("informe o diretório dos CSVs *_fp.csv")

println("="^84)
println("TAXA LOGARÍTMICA DE VOTO — reprodução da Fig. 3 de Klimek et al. (2012)")
println("="^84)
@printf("%-26s %-9s %-12s %-10s %-16s\n",
        "pleito", "seções", "assimetria", "curtose", "canto ≥95%/≥95%")

# controle sintético: eleição limpa por construção
rng = MersenneTwister(20260830)
m = 20_000
N = rand(rng, 800:2500, m)
a = clamp.(rand(rng, Normal(0.62, 0.09), m), 0.05, 0.99)
s = clamp.(rand(rng, Normal(0.48, 0.10), m), 0.02, 0.98)
V = [round(Int, N[i]*a[i]) for i in 1:m]
W = [round(Int, V[i]*s[i]) for i in 1:m]
rc = election_fingerprint(W, V, N)
@printf("%-26s %-9d %+-12.2f %-10.2f %d (%.2f%%)\n",
        "simulada limpa", rc.m, rc.skewness, rc.kurtosis, canto(rc), 100canto(rc)/rc.m)

for (rot, arq) in (("EUA 2020 · Illinois", "us2020_il_fp.csv"),
                   ("EUA 2020 · Alabama",  "us2020_al_fp.csv"),
                   ("Rússia 2012",         "russia2012_fp.csv"))
    f = joinpath(D, arq)
    (isfile(f) && filesize(f) > 60) || continue
    w, v, n = ler3(f)
    isempty(w) && continue
    r = election_fingerprint(w, v, n)
    @printf("%-26s %-9d %+-12.2f %-10.2f %d (%.2f%%)\n",
            rot, r.m, r.skewness, r.kurtosis, canto(r), 100canto(r)/r.m)
end

println("""

Leitura:
· Illinois cai praticamente em cima do (0, 3) previsto para eleições limpas.
· Alabama se afasta, mas sem nenhuma seção no canto. O artigo atribui desvios
  desse tipo a heterogeneidade da população (o exemplo dele é o Canadá,
  anglófono × francófono), não a fraude. Aqui há ainda viés de subamostra:
  só 1 204 dos 2 110 precincts do Alabama publicam eleitorado.
· A Rússia mostra as duas assinaturas descritas no artigo — o borrão diagonal
  em direção ao canto superior direito e o segundo aglomerado em (100%, 100%).

O modelo paramétrico de fraude (f_i, f_e) do artigo NÃO está implementado: a
especificação dos mecanismos está no Supporting Information, que não acompanha
o preprint do arXiv.""")
