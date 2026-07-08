# ── último dígito (Beber & Scacco 2012) ──────────────────────────────

struct LastDigitResult
    position::Symbol
    counts::Vector{Int}
    observed::Vector{Float64}
    n::Int
    n_excluded::Int
    chi2::Float64
    df::Int
    pvalue::Float64
    freq_0_5::Float64
end

"""
    last_digit_test(x::AbstractVector{<:Integer};
                    position::Symbol = :last,
                    min_value::Integer = position === :last ? 10 : 100)
        -> LastDigitResult

Testa a uniformidade do último (`position = :last`) ou penúltimo
(`position = :penultimate`) dígito de contagens de votos, via χ² contra
a distribuição uniforme em 0:9 (Beber & Scacco 2012).

Em dados limpos, o último dígito de contagens grandes é ~uniforme;
números fabricados por humanos superusam 0 e 5 e subusam repetições.
`min_value` exclui contagens pequenas, cujo último dígito não é uniforme
mesmo em dados limpos (default: 10 para o último, 100 para o penúltimo).

O campo `freq_0_5` traz a frequência conjunta dos dígitos 0 e 5
(esperado ≈ 0.20).
"""
function last_digit_test(x::AbstractVector{<:Integer};
                         position::Symbol = :last,
                         min_value::Integer = position === :last ? 10 : 100)
    position in (:last, :penultimate) ||
        throw(ArgumentError("position deve ser :last ou :penultimate"))
    getd = position === :last ? last_digit : penultimate_digit

    v = [Int(xi) for xi in x if xi ≥ min_value]
    n = length(v)
    n_excluded = length(x) - n
    n > 0 || throw(ArgumentError("nenhuma observação válida (x ≥ $min_value)"))
    n < 100 && @warn "amostra pequena para teste de dígito (n = $n)"

    counts = zeros(Int, 10)
    @inbounds for xi in v
        counts[getd(xi) + 1] += 1
    end

    observed = counts ./ n
    E = n / 10
    chi2 = sum((counts .- E) .^ 2 ./ E)
    df = 9
    pvalue = ccdf(Chisq(df), chi2)
    freq_0_5 = observed[1] + observed[6]

    return LastDigitResult(position, counts, observed, n, n_excluded,
                           chi2, df, pvalue, freq_0_5)
end

function Base.show(io::IO, ::MIME"text/plain", r::LastDigitResult)
    nome = r.position === :last ? "último" : "penúltimo"
    printstyled(io, "Teste de uniformidade — $nome dígito\n"; bold = true)
    println(io, "n = $(r.n) (excluídos: $(r.n_excluded))")
    pmax = maximum(r.observed)
    for d in 0:9
        print(io, "  ", d, " │", _bar(r.observed[d + 1], pmax), " ")
        println(io, _fmt_pct(r.observed[d + 1]), "  (esp.  10.0%)")
    end
    @printf(io, "χ²(%d) = %.2f   p = %s   freq{0,5} = %s (esp. 20.0%%)\n",
            r.df, r.chi2, _fmt_p(r.pvalue), _fmt_pct(r.freq_0_5))
    print(io, "H₀ (uniforme): ")
    printstyled(io, r.pvalue < 0.05 ? "rejeitada" : "não rejeitada";
                color = _sig_color(r.pvalue), bold = true)
    println(io)
end
