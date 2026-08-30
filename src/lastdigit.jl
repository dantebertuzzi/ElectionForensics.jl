# ── último dígito (Beber & Scacco 2012) ──────────────────────────────

"""
    LastDigitResult

Resultado de [`last_digit_test`](@ref).

Campos: `position` (`:last` ou `:penultimate`), `counts`, `observed`, `n`,
`n_excluded`, `chi2`, `df`, `pvalue`, `freq_0_5` (frequência conjunta dos
dígitos 0 e 5, esperada ≈ 0,20) e `null_valid` — `false` quando o null uniforme
**não** é defensável para os dados fornecidos (veja [`last_digit_test`](@ref)).
"""
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
    null_valid::Bool
end

"""
    DigitTestConfig

Parâmetros do teste de uniformidade de dígitos.

- `position`: `:last` (recomendado) ou `:penultimate`.
- `min_value`: contagens abaixo deste valor são excluídas. O default de 10 para
  `:last` é o de Beber & Scacco; contagens muito pequenas (`sd < ~4`) quebram a
  uniformidade — veja `validation/02_calibracao_ultimodigito.jl`.
- `min_n`: abaixo disso é emitido aviso de amostra pequena.
- `min_expected`: regra de Cochran para a validade da aproximação χ².
- `warn`: emite avisos de aplicabilidade.
"""
Base.@kwdef struct DigitTestConfig
    position::Symbol = :last
    min_value::Int = position === :last ? 10 : 100
    min_n::Int = 100
    min_expected::Float64 = 5.0
    warn::Bool = true
end

# ordens de grandeza cobertas pelas contagens; o null uniforme do penúltimo
# dígito exige densidade suave numa escala de 100 — na prática, ≥ 2 décadas.
_orders_of_magnitude(v) =
    log10(quantile(Float64.(v), 0.95) / max(quantile(Float64.(v), 0.05), 1.0))

"""
    last_digit_test(x::AbstractVector{<:Integer}; position = :last, kwargs...)
        -> LastDigitResult
    last_digit_test(x, cfg::DigitTestConfig) -> LastDigitResult

Testa a uniformidade do último (`position = :last`) ou penúltimo
(`position = :penultimate`) dígito de contagens de votos, via χ² contra a
distribuição uniforme em 0:9 (Beber & Scacco 2012).

O campo `freq_0_5` traz a frequência conjunta dos dígitos 0 e 5 (esperado ≈ 0,20):
números fabricados por humanos superusam 0 e 5.

!!! warning "O null uniforme do PENÚLTIMO dígito raramente é válido"
    A uniformidade do último dígito exige apenas que a densidade das contagens
    seja suave numa escala de 10 — condição satisfeita em dados eleitorais
    reais (erro tipo I ≈ α para seções com ≥ 60 eleitores; verificado em
    `validation/02_calibracao_ultimodigito.jl`).

    O **penúltimo** dígito exige suavidade numa escala de 100, o que quase nunca
    vale: em eleições limpas simuladas com seções de 150–900 eleitores o teste
    rejeita H₀ em 50 % dos casos, e em 100 % com seções de 100–200
    (`validation/07_penultimo_digito.jl`). Sob dados limpos a distribuição do
    penúltimo dígito é **decrescente** (10,9 % → 8,9 %), não uniforme. Não há
    guarda computável a partir dos dados que salve esse null — testamos um null
    reamostrado e ele reduz mas não elimina a inflação. Trate `:penultimate`
    como diagnóstico exploratório, nunca como teste; `null_valid` sinaliza
    quando a condição mínima (≥ 2 décadas de amplitude) não é atendida.
"""
function last_digit_test(x::AbstractVector{<:Integer}, cfg::DigitTestConfig)
    cfg.position in (:last, :penultimate) ||
        throw(ArgumentError("position deve ser :last ou :penultimate"))
    # Abaixo do mínimo estrutural o dígito não existe (penúltimo de x < 100) ou
    # não pode ser uniforme nem em dados limpos.
    minimo = cfg.position === :last ? 10 : 100
    cfg.min_value ≥ minimo ||
        throw(ArgumentError("min_value = $(cfg.min_value) é menor que o mínimo \
            estrutural de :$(cfg.position) ($minimo)"))
    getd = cfg.position === :last ? last_digit : penultimate_digit

    v = [Int(xi) for xi in x if xi ≥ cfg.min_value]
    n = length(v)
    n_excluded = length(x) - n
    n > 0 || throw(ArgumentError("nenhuma observação válida (x ≥ $(cfg.min_value))"))

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

    oom = _orders_of_magnitude(v)
    null_valid = cfg.position === :last ? true : oom ≥ 2.0

    if cfg.warn
        n < cfg.min_n && @warn "amostra pequena para teste de dígito (n = $n)"
        # A uniformidade do último dígito exige sd(contagens) >> 10; com seções
        # muito pequenas o erro tipo I sobe a 0,30 (validation/02).
        sdv = length(v) > 1 ? std(v) : 0.0
        if cfg.position === :last && sdv < 4 * 10
            @warn "contagens pouco dispersas (sd = $(round(sdv, digits = 1))) \
                para um null uniforme em base 10; considere `min_value` maior \
                ou trate o p-valor como descritivo"
        end
        E < cfg.min_expected && @warn "aproximação χ² inválida: casela esperada \
            = $(round(E, digits = 2)) < $(cfg.min_expected)"
        if cfg.position === :penultimate && !null_valid
            @warn "null uniforme do penúltimo dígito NÃO é válido para estes \
                dados: contagens cobrem $(round(oom, digits = 2)) décadas (< 2). \
                Em eleições limpas com essa amplitude o teste rejeita H₀ em até \
                100 % dos casos. Trate o p-valor como descritivo."
        end
    end

    return LastDigitResult(cfg.position, counts, observed, n, n_excluded,
                           chi2, df, pvalue, freq_0_5, null_valid)
end

last_digit_test(x::AbstractVector{<:Integer}; kwargs...) =
    last_digit_test(x, DigitTestConfig(; kwargs...))

last_digit_test(x::AbstractVector; skipmissing::Bool = false, kwargs...) =
    last_digit_test(_coerce_counts(x; skipmissing = skipmissing, name = "x"),
                    DigitTestConfig(; kwargs...))

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
    if !r.null_valid
        printstyled(io, "⚠ null uniforme não é válido para estes dados — \
                        p-valor é descritivo, não um teste\n"; color = :red)
    end
    print(io, "H₀ (uniforme): ")
    printstyled(io, r.pvalue < 0.05 ? "rejeitada" : "não rejeitada";
                color = _sig_color(r.pvalue), bold = true)
    println(io)
end
