# ── extração de dígitos ──────────────────────────────────────────────

"""
    first_digit(x::Integer) -> Int

Primeiro dígito significativo de `x`. Requer `x ≥ 1`.
"""
function first_digit(x::Integer)
    x ≥ 1 || throw(ArgumentError("first_digit requer x ≥ 1 (recebido: $x)"))
    while x ≥ 10
        x ÷= 10
    end
    return Int(x)
end

"""
    second_digit(x::Integer) -> Int

Segundo dígito significativo de `x`. Requer `x ≥ 10`.
"""
function second_digit(x::Integer)
    x ≥ 10 || throw(ArgumentError("second_digit requer x ≥ 10 (recebido: $x)"))
    while x ≥ 100
        x ÷= 10
    end
    return Int(x % 10)
end

"""
    last_digit(x::Integer) -> Int

Último dígito (unidades) de `|x|`.
"""
last_digit(x::Integer) = Int(abs(x) % 10)

"""
    penultimate_digit(x::Integer) -> Int

Penúltimo dígito (dezenas) de `|x|`. Requer `|x| ≥ 10`.
"""
function penultimate_digit(x::Integer)
    a = abs(x)
    a ≥ 10 || throw(ArgumentError("penultimate_digit requer |x| ≥ 10 (recebido: $x)"))
    return Int((a ÷ 10) % 10)
end

# ── helpers internos de exibição ─────────────────────────────────────

const _BAR_WIDTH = 24

function _bar(p::Real, pmax::Real; width::Int = _BAR_WIDTH)
    pmax > 0 || return " "^width
    n = clamp(round(Int, width * p / pmax), 0, width)
    return "█"^n * "░"^(width - n)
end

_fmt_pct(p) = @sprintf("%5.1f%%", 100p)
_fmt_p(p) = p < 0.001 ? @sprintf("%.2e", p) : @sprintf("%.4f", p)

function _sig_color(p)
    p < 0.01  && return :red
    p < 0.05  && return :yellow
    return :green
end

# Piso de banda usado quando a dispersão amostral é nula ou não finita — por
# exemplo, todas as seções com o mesmo percentual. Não é uma escolha
# estatística: é o menor jitter que ainda produz uma nula não degenerada.
# `_silverman` sinaliza esse caso devolvendo `degenerate = true`.
const _H_FLOOR = 0.01

"""
    _silverman(v) -> (h, degenerate)

Banda de Silverman `0.9 · min(σ, IQR/1.34) · m^(-1/5)`. Quando a dispersão é
nula ou não finita, devolve `(_H_FLOOR, true)` — o chamador decide se avisa.
"""
function _silverman(v::AbstractVector{<:Real})
    m = length(v)
    m > 1 || return (_H_FLOOR, true)
    s = std(v)
    iqrv = quantile(v, 0.75) - quantile(v, 0.25)
    spread = min(s, iqrv / 1.34)
    spread ≤ 0 && (spread = max(s, iqrv / 1.34))
    h = 0.9 * spread * m^(-1 / 5)
    return (isfinite(h) && h > 0) ? (h, false) : (_H_FLOOR, true)
end

# ajuste de Benjamini–Hochberg (FDR) para p-valores múltiplos
function _bh_adjust(p::AbstractVector{<:Real})
    n = length(p)
    n == 0 && return Float64[]
    ord = sortperm(p)
    q = Vector{Float64}(undef, n)
    running = 1.0
    @inbounds for k in n:-1:1
        running = min(running, n * p[ord[k]] / k)
        q[ord[k]] = running
    end
    return q
end


# ── coerção de entrada ───────────────────────────────────────────────
#
# CSVs eleitorais chegam como Float64 ou com `missing`. Coagimos de forma
# explícita e ruidosa: nada é descartado em silêncio numa ferramenta forense.

_coerce_counts(x::AbstractVector{<:Integer}; kwargs...) = x

function _coerce_counts(x::AbstractVector; skipmissing::Bool = false,
                        name::AbstractString = "x")
    nmiss = count(ismissing, x)
    if nmiss > 0 && !skipmissing
        throw(ArgumentError("$name contém $nmiss valores `missing`; passe \
            `skipmissing = true` para descartá-los explicitamente"))
    end
    out = Vector{Int}(undef, 0)
    sizehint!(out, length(x) - nmiss)
    @inbounds for xi in x
        ismissing(xi) && continue
        xi isa Real || throw(ArgumentError("$name contém elemento não numérico: $xi"))
        isfinite(xi) || throw(ArgumentError("$name contém `$xi`; contagens devem ser finitas"))
        isinteger(xi) || throw(ArgumentError("$name contém o valor não inteiro $xi; \
            contagens de votos devem ser inteiras"))
        push!(out, Int(xi))
    end
    return out
end
