# ── Lei de Benford (1BL / 2BL) ───────────────────────────────────────

"""
    benford_expected(digit::Int) -> Vector{Float64}

Probabilidades esperadas sob a lei de Benford para o 1º dígito
(d = 1:9) ou 2º dígito (d = 0:9).
"""
function benford_expected(digit::Int)
    if digit == 1
        return [log10(1 + 1 / d) for d in 1:9]
    elseif digit == 2
        return [sum(log10(1 + 1 / (10k + d)) for k in 1:9) for d in 0:9]
    else
        throw(ArgumentError("digit deve ser 1 ou 2"))
    end
end

# limiares de MAD (Nigrini 2012)
const _NIGRINI = Dict(
    1 => (0.006, 0.012, 0.015),
    2 => (0.008, 0.010, 0.012),
)

function _nigrini_conformity(mad::Float64, digit::Int)
    t1, t2, t3 = _NIGRINI[digit]
    mad < t1 && return :conforme
    mad < t2 && return :aceitavel
    mad < t3 && return :marginal
    return :nao_conforme
end

struct BenfordResult
    digit::Int
    labels::Vector{Int}
    counts::Vector{Int}
    observed::Vector{Float64}
    expected::Vector{Float64}
    n::Int
    n_excluded::Int
    chi2::Float64
    df::Int
    pvalue::Float64
    mad::Float64
    conformity::Symbol
end

"""
    benford_test(x::AbstractVector{<:Integer}; digit::Int = 2) -> BenfordResult

Testa a distribuição do `digit`-ésimo dígito significativo de `x` contra
a lei de Benford, via χ² de aderência e MAD (com limiares de conformidade
de Nigrini).

- `digit = 1` (1BL) exige `x ≥ 1`; `digit = 2` (2BL, Mebane 2008) exige
  `x ≥ 10`. Valores abaixo do mínimo são excluídos (`n_excluded`).
- Para dados eleitorais, prefira `digit = 2`: a aderência do 1º dígito em
  contagens de votos é contestada (Deckert, Myagkov & Ordeshook 2011),
  pois depende da distribuição de tamanhos das seções.

Retorna um [`BenfordResult`](@ref) com frequências observadas/esperadas,
χ², p-valor, MAD e classificação de conformidade
(`:conforme`, `:aceitavel`, `:marginal`, `:nao_conforme`).
"""
function benford_test(x::AbstractVector{<:Integer}; digit::Int = 2)
    haskey(_NIGRINI, digit) || throw(ArgumentError("digit deve ser 1 ou 2"))
    minv = digit == 1 ? 1 : 10
    v = [Int(xi) for xi in x if xi ≥ minv]
    n = length(v)
    n_excluded = length(x) - n
    n > 0 || throw(ArgumentError("nenhuma observação válida (x ≥ $minv)"))
    n < 100 && @warn "amostra pequena para teste de Benford (n = $n)"

    labels = digit == 1 ? collect(1:9) : collect(0:9)
    getd = digit == 1 ? first_digit : second_digit
    counts = zeros(Int, length(labels))
    offset = first(labels)
    @inbounds for xi in v
        counts[getd(xi) - offset + 1] += 1
    end

    expected = benford_expected(digit)
    observed = counts ./ n
    E = n .* expected
    chi2 = sum((counts .- E) .^ 2 ./ E)
    df = length(labels) - 1
    pvalue = ccdf(Chisq(df), chi2)
    mad = mean(abs.(observed .- expected))

    return BenfordResult(digit, labels, counts, observed, expected,
                         n, n_excluded, chi2, df, pvalue, mad,
                         _nigrini_conformity(mad, digit))
end

function Base.show(io::IO, ::MIME"text/plain", r::BenfordResult)
    printstyled(io, "Teste de Benford — $(r.digit)º dígito ($(r.digit)BL)\n";
                bold = true)
    println(io, "n = $(r.n) (excluídos: $(r.n_excluded))")
    pmax = max(maximum(r.observed), maximum(r.expected))
    for (i, d) in enumerate(r.labels)
        print(io, "  ", d, " │", _bar(r.observed[i], pmax), " ")
        print(io, _fmt_pct(r.observed[i]))
        println(io, "  (esp. ", _fmt_pct(r.expected[i]), ")")
    end
    @printf(io, "χ²(%d) = %.2f   p = %s   MAD = %.4f\n",
            r.df, r.chi2, _fmt_p(r.pvalue), r.mad)
    print(io, "conformidade (Nigrini): ")
    printstyled(io, String(r.conformity);
                color = r.conformity in (:conforme, :aceitavel) ? :green :
                        r.conformity === :marginal ? :yellow : :red,
                bold = true)
    println(io)
end
