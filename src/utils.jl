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

# largura de banda de Silverman com salvaguardas
function _silverman(v::AbstractVector{<:Real})
    m = length(v)
    m > 1 || return 0.01
    s = std(v)
    iqrv = quantile(v, 0.75) - quantile(v, 0.25)
    spread = min(s, iqrv / 1.34)
    spread ≤ 0 && (spread = max(s, iqrv / 1.34))
    h = 0.9 * spread * m^(-1 / 5)
    return (isfinite(h) && h > 0) ? h : 0.01
end
