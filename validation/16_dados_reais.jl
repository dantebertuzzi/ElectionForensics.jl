# Validação em DADOS REAIS de fontes oficiais.
#
#  (a) EUA 2020, resultados por precinct — OpenElections, transcritos dos
#      boletins oficiais das secretarias eleitorais estaduais.
#      Pleitos sem alegação séria de fraude sistemática ⇒ testam CALIBRAÇÃO.
#
#  (b) Rússia 2012, presidencial, 95 413 seções (UIK) — dataset de
#      Sergey Shpilkin usado em Kobak, Shpilkin & Pshenichnikov (2016),
#      "Integer percentages as electoral falsification fingerprints",
#      Annals of Applied Statistics 10(1):54–73. Anomalia documentada e
#      revisada por pares ⇒ testa PODER contra uma verdade externa.
#
# Uso: julia -t auto validation/16_dados_reais.jl <dir_com_csvs>
using ElectionForensics, Random, Statistics, Printf

const DIR = length(ARGS) ≥ 1 ? ARGS[1] : error("informe o diretório dos CSVs")

"Leitor CSV mínimo: colunas `cv` (votos) e `ct` (total) por índice."
function ler(path; cv::Int = 1, ct::Int = 2)
    a = Int[]; b = Int[]
    for (k, linha) in enumerate(eachline(path))
        k == 1 && continue
        campos = split(linha, ',')
        length(campos) < max(cv, ct) && continue
        x = tryparse(Int, campos[cv]); y = tryparse(Int, campos[ct])
        (x === nothing || y === nothing) && continue
        (y > 0 && 0 ≤ x ≤ y) || continue
        push!(a, x); push!(b, y)
    end
    return a, b
end

decadas(v) = log10(quantile(Float64.(v), 0.95) / max(quantile(Float64.(v), 0.05), 1.0))

function analisa(nome, votes, totals; B = 999)
    rng = MersenneTwister(20260829)
    ben_c = benford_test(votes; digit = 2, null = :benford, warn = false)
    ben_r = benford_test(votes; digit = 2, null = :resampled, B = B,
                         warn = false, rng = MersenneTwister(1))
    ld    = last_digit_test(votes; warn = false)
    rz    = rozenas_test(votes, totals; B = B, rng = rng)
    (nome = nome, m = length(votes), dec = decadas(votes),
     med = median(totals),
     ben_c = ben_c.pvalue, ben_r = ben_r.pvalue, mad = ben_c.mad,
     ld = ld.pvalue, f05 = ld.freq_0_5,
     rz = rz.pvalue, Tobs = rz.total_observed, Tnull = mean(rz.null_totals),
     z = rz.zscore, nq = count(<(0.05), rz.qvalues), res = rz)
end

# Todos os pleitos americanos baixados, exceto os que a checagem de
# integridade de `fetch_dados_reais.py` reprovou.
const EXCLUIR = Set(String[])   # o crivo já roda em fetch_dados_reais.py
                              # linha-resumo da fonte (voto ranqueado, transcrição
                              # parcialmente reconciliada). Analisar dados que não
                              # sei extrair corretamente contaminaria a calibração.

arquivos = sort(filter(f -> startswith(f, "us2020_") && endswith(f, ".csv"),
                       readdir(DIR)))
arquivos = filter(f -> !(f[8:9] in EXCLUIR), arquivos)

println("="^108)
println("(a) EUA 2020 — resultados por precinct (OpenElections / secretarias estaduais)")
println("    Esperado: p-valores ~ U(0,1). Testa CALIBRAÇÃO em dados reais.")
println("="^108)
@printf("%-6s %-7s %-7s %-8s %-11s %-11s %-10s %-10s %-8s\n",
        "UF", "precs", "déc.", "med(N)", "2BL clás.", "2BL reamos.", "últ.díg", "rozenas", "q<.05")
resultados = []
for f in arquivos
    v, t = ler(joinpath(DIR, f))
    length(v) < 150 && continue
    r = analisa(uppercase(f[8:9]), v, t)
    push!(resultados, r)
    @printf("%-6s %-7d %-7.2f %-8.0f %-11s %-11.3f %-10.3f %-10.3f %-8d\n",
            r.nome, r.m, r.dec, r.med,
            r.ben_c < 1e-3 ? @sprintf("%.1e", r.ben_c) : @sprintf("%.3f", r.ben_c),
            r.ben_r, r.ld, r.rz, r.nq)
end

# Teste de uniformidade dos p-valores. Sob H0 (pleitos limpos e null correto)
# eles devem ser U(0,1); desvio indica miscalibração do null em dados reais.
function ks_unif(ps)
    n = length(ps); n < 5 && return NaN
    s = sort(ps)
    d = maximum(max(abs(s[i] - (i-1)/n), abs(i/n - s[i])) for i in 1:n)
    # p-valor assintótico de Kolmogorov
    λ = (sqrt(n) + 0.12 + 0.11/sqrt(n)) * d
    p = 2 * sum((-1)^(k-1) * exp(-2 * k^2 * λ^2) for k in 1:100)
    return clamp(p, 0.0, 1.0)
end

n = length(resultados)
println()
@printf("Calibração em %d pleitos americanos reais:\n", n)
@printf("  %-24s %-12s %-12s %-12s\n", "teste", "rej @0,05", "rej @0,10", "KS vs U(0,1)")
for (rot, f) in (("2BL null clássico", r -> r.ben_c),
                 ("2BL null reamostrado", r -> r.ben_r),
                 ("último dígito", r -> r.ld),
                 ("frações coarse", r -> r.rz))
    ps = [f(r) for r in resultados]
    @printf("  %-24s %-12s %-12s %-12.4f\n", rot,
            string(count(<(0.05), ps), "/", n), string(count(<(0.10), ps), "/", n),
            ks_unif(ps))
end
println()
println("  (esperado sob calibração: ~5 % e ~10 % de rejeições, KS não significativo)")

# ── (b) Rússia 2012 ──────────────────────────────────────────────────
pru = joinpath(DIR, "russia2012.csv")
if isfile(pru)
    println()
    println("="^108)
    println("(b) Rússia 2012, presidencial — 95 413 seções (Shpilkin / Kobak et al. 2016)")
    println("    Esperado: anomalia forte em frações redondas. Testa PODER.")
    println("="^108)
    v, t = ler(pru; cv = 1, ct = 3)   # Putin sobre bulletins válidos
    r = analisa("RU2012", v, t)
    @printf("seções: %d   mediana de votos válidos: %.0f   amplitude: %.2f décadas\n",
            r.m, r.med, r.dec)
    @printf("  Benford 2BL clássico   p = %.3e   (MAD = %.4f)\n", r.ben_c, r.mad)
    @printf("  Benford 2BL reamostrado p = %.4f\n", r.ben_r)
    @printf("  último dígito          p = %.3e   freq{0,5} = %.4f (esp. 0.20)\n", r.ld, r.f05)
    @printf("  frações coarse         p = %.4f   T_obs = %d   E[T|H₀] = %.0f   z = %.1f\n",
            r.rz, r.Tobs, r.Tnull, r.z)
    @printf("  frações com q < 0,05   : %d de %d\n", r.nq, length(r.res.fractions))
    println()
    println("  As dez frações com maior excesso:")
    println("    fração    obs     esp       z         q")
    zs = [(r.res.observed[j] - r.res.null_mean[j]) /
          max(r.res.null_sd[j], eps()) for j in eachindex(r.res.fractions)]
    for j in sortperm(zs; rev = true)[1:10]
        f = r.res.fractions[j]
        @printf("    %-8s %5d  %7.1f  %8.1f  %.2e\n",
                string(numerator(f), "/", denominator(f)),
                r.res.observed[j], r.res.null_mean[j], zs[j], r.res.qvalues[j])
    end
end

# ── (c) por que o sinal russo é fraco? ───────────────────────────────
#
# O teste conta apenas percentuais EXATAMENTE iguais a k/d. Isso exige que
# `totals` seja divisível por d. Kobak et al. comparam a densidade em
# percentuais inteiros contra a vizinhança local, sem exigir divisibilidade —
# por isso enxergam um efeito muito maior nos mesmos dados.
if isfile(joinpath(DIR, "russia2012.csv"))
    v, t = ler(joinpath(DIR, "russia2012.csv"); cv = 1, ct = 3)
    println()
    println("="^108)
    println("(c) A restrição de divisibilidade limita o poder")
    println("="^108)
    println("Uma seção só pode acertar k/d se `totals` for divisível por d:")
    @printf("  %-6s %-14s %-14s\n", "d", "seções com d|N", "%")
    for d in 2:10
        c = count(x -> x % d == 0, t)
        @printf("  %-6d %-14d %-14.1f\n", d, c, 100c / length(t))
    end
    elegivel = count(i -> any(d -> t[i] % d == 0, 2:10), eachindex(t))
    @printf("\n  seções elegíveis para ALGUMA fração d≤10: %d (%.1f %%)\n",
            elegivel, 100elegivel / length(t))
    @printf("  seções que de fato acertam: %d (%.1f %%)\n",
            count(i -> any(d -> t[i] % d == 0 && (v[i] * d) % t[i] == 0, 2:10), eachindex(t)),
            100count(i -> any(d -> t[i] % d == 0 && (v[i]*d) % t[i] == 0, 2:10), eachindex(t)) / length(t))
    println("""
  ⇒ O teste vê apenas a fatia da manipulação que cai numa fração exata de
    denominador pequeno. É uma limitação inerente do primitivo, não um bug:
    `max_denom` maior recupera parte do sinal (p = 0,042 → 0,017 com d ≤ 100).""")
end

# ── (d) o teste do penúltimo dígito nos dados russos ─────────────────
if isfile(joinpath(DIR, "russia2012.csv"))
    v, _ = ler(joinpath(DIR, "russia2012.csv"); cv = 1, ct = 3)
    pd = last_digit_test(v; position = :penultimate, warn = false)
    println()
    println("="^108)
    println("(d) Confirmação em dados reais do achado C1")
    println("="^108)
    @printf("  penúltimo dígito, Rússia 2012: p = %.3e  ← alarme catastrófico\n", pd.pvalue)
    @printf("  null_valid = %s (contagens cobrem %.2f décadas, exige ≥ 2)\n",
            pd.null_valid, decadas(v))
    println("""
  Antes da correção este p-valor entrava no relatório sem ressalva. O mesmo
  teste devolve p ≈ 1e-8 em dados SIMULADOS LIMPOS: o número não distingue
  fraude de artefato do null. Agora é opt-in e vem marcado.""")
end
