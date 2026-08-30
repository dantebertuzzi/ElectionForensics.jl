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

# limiares de MAD (Nigrini 2012, Tab. 5.1) — heurísticos, não distribucionais
const _NIGRINI = Dict(
    1 => (0.006, 0.012, 0.015),
    2 => (0.008, 0.010, 0.012),
)

# Os intervalos de Nigrini (2012, Tab. 5.1) são fechados no limite superior:
# 0,000–0,006 conforme; 0,006–0,012 aceitável; 0,012–0,015 marginal.
function _nigrini_conformity(mad::Float64, digit::Int)
    t1, t2, t3 = _NIGRINI[digit]
    mad ≤ t1 && return :conforme
    mad ≤ t2 && return :aceitavel
    mad ≤ t3 && return :marginal
    return :nao_conforme
end

"""
    BenfordConfig

Parâmetros do teste de Benford.

- `digit`: 1 (1BL) ou 2 (2BL, Mebane 2008).
- `null`: distribuição de referência do p-valor.
  - `:resampled` (default) — null **reamostrado**: as contagens são reamostradas
    com jitter gaussiano de banda `log_bandwidth` em `log10`, larga o bastante
    para apagar a estrutura do dígito mas estreita para preservar o formato da
    distribuição de contagens. Responde "o dígito é anômalo **dado** o formato
    empírico das contagens?".
  - `:benford` — χ² assintótico contra a lei de Benford pura. **Não é calibrado
    para contagens eleitorais**: em seções de tamanho homogêneo a taxa de erro
    tipo I chega a 0,72 (α = 0,05). Use apenas como estatística descritiva.
- `B`: réplicas do null reamostrado.
- `log_bandwidth`: banda do kernel, em décadas.
- `min_n`, `min_expected`: guardas de validade (tamanho amostral e regra de
  Cochran para a aproximação χ²).
- `warn`: emite avisos de aplicabilidade.
"""
Base.@kwdef struct BenfordConfig
    digit::Int = 2
    null::Symbol = :resampled
    B::Int = 999
    log_bandwidth::Float64 = 0.05
    min_n::Int = 100
    min_expected::Float64 = 5.0
    warn::Bool = true
end

"""
    BenfordResult

Resultado de [`benford_test`](@ref).

Campos: `digit`, `labels` (dígitos testados), `counts`, `observed` e `expected`
(frequências relativas), `n`, `n_excluded`, `chi2`, `df`, `pvalue` (segundo
`null`), `pvalue_asymptotic` (χ² contra Benford puro, sempre presente), `null`,
`mad` e `conformity` (classificação de Nigrini).

`pvalue` só é interpretável como p-valor de teste quando `null === :resampled`;
veja [`BenfordConfig`](@ref).
"""
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
    pvalue_asymptotic::Float64
    null::Symbol
    mad::Float64
    conformity::Symbol
end

_first_digit_of(x) = (while x ≥ 10; x ÷= 10; end; Int(x))
_second_digit_of(x) = (while x ≥ 100; x ÷= 10; end; Int(x % 10))

# `getd` como valor de tipo (Val) em vez de ternário: sem Union de tipos de
# função no laço quente.
function _digit_counts!(counts::Vector{Int}, v::AbstractVector{Int},
                        getd::F, offset::Int) where {F}
    fill!(counts, 0)
    @inbounds for xi in v
        counts[getd(xi) - offset + 1] += 1
    end
    return counts
end

_digit_counts(v, digit::Int, ncat::Int, offset::Int) =
    digit == 1 ? _digit_counts!(zeros(Int, ncat), v, _first_digit_of, offset) :
                 _digit_counts!(zeros(Int, ncat), v, _second_digit_of, offset)

_chi2(counts, E) = sum((counts .- E) .^ 2 ./ E)

# null reamostrado: jitter gaussiano em log10 apaga a estrutura do dígito
# preservando o formato macro da distribuição de contagens.
function _resampled_pvalue(v::AbstractVector{Int}, digit::Int, expected::Vector{Float64},
                           chi_obs::Float64, cfg::BenfordConfig, rng::AbstractRNG)
    n = length(v)
    minv = digit == 1 ? 1 : 10
    lv = log10.(Float64.(v))
    ncat, offset = length(expected), (digit == 1 ? 1 : 0)
    E = n .* expected
    ge = 0
    buf = Vector{Int}(undef, n)
    @inbounds for _ in 1:cfg.B
        k = 0
        while k < n
            y = round(Int, 10^(lv[rand(rng, 1:n)] + cfg.log_bandwidth * randn(rng)))
            if y ≥ minv
                k += 1
                buf[k] = y
            end
        end
        _chi2(_digit_counts(buf, digit, ncat, offset), E) ≥ chi_obs && (ge += 1)
    end
    return (1 + ge) / (cfg.B + 1)
end

"""
    benford_test(x::AbstractVector{<:Integer}; digit = 2, null = :resampled,
                 rng = Random.default_rng(), kwargs...) -> BenfordResult
    benford_test(x, cfg::BenfordConfig; rng = Random.default_rng()) -> BenfordResult

Testa a distribuição do `digit`-ésimo dígito significativo de `x`.

- `digit = 1` (1BL) exige `x ≥ 1`; `digit = 2` (2BL) exige `x ≥ 10`. Valores
  abaixo do mínimo são excluídos (`n_excluded`).

!!! warning "Contagens eleitorais não seguem a lei de Benford"
    Sob eleições limpas simuladas, o χ² contra a lei de Benford rejeita H₀ em
    até 100 % dos casos (1BL) e 72 % (2BL) quando as seções têm tamanho
    homogêneo — o caso brasileiro. A conformidade só vale quando os tamanhos de
    seção são aproximadamente log-normais com dispersão grande
    (Deckert, Myagkov & Ordeshook 2011). Por isso o default é `null = :resampled`,
    que calibra o p-valor contra a própria distribuição empírica de contagens.
    `null = :benford` reproduz o teste clássico e deve ser lido como estatística
    **descritiva** (MAD/conformidade), não como teste de hipótese.

Retorna um [`BenfordResult`](@ref). Veja [`BenfordConfig`](@ref) para os
parâmetros.
"""
function benford_test(x::AbstractVector{<:Integer}, cfg::BenfordConfig;
                      rng::AbstractRNG = Random.default_rng())
    haskey(_NIGRINI, cfg.digit) || throw(ArgumentError("digit deve ser 1 ou 2"))
    cfg.null in (:resampled, :benford) ||
        throw(ArgumentError("null deve ser :resampled ou :benford"))
    cfg.B ≥ 99 || throw(ArgumentError("B deve ser ≥ 99"))
    cfg.log_bandwidth > 0 || throw(ArgumentError("log_bandwidth deve ser > 0"))

    minv = cfg.digit == 1 ? 1 : 10
    v = [Int(xi) for xi in x if xi ≥ minv]
    n = length(v)
    n_excluded = length(x) - n
    n > 0 || throw(ArgumentError("nenhuma observação válida (x ≥ $minv)"))

    expected = benford_expected(cfg.digit)
    labels = cfg.digit == 1 ? collect(1:9) : collect(0:9)
    offset = first(labels)
    counts = _digit_counts(v, cfg.digit, length(labels), offset)
    observed = counts ./ n
    E = n .* expected
    chi2 = _chi2(counts, E)
    df = length(labels) - 1
    p_asym = ccdf(Chisq(df), chi2)
    mad = mean(abs.(observed .- expected))

    if cfg.warn
        n < cfg.min_n && @warn "amostra pequena para teste de Benford (n = $n)"
        # regra de Cochran: a aproximação χ² exige todas as caselas esperadas ≥ 5
        minE = minimum(E)
        minE < cfg.min_expected && @warn "aproximação χ² inválida: menor casela \
            esperada = $(round(minE, digits = 2)) < $(cfg.min_expected) (n = $n)"
        if cfg.null === :benford
            @warn "null = :benford não é calibrado para contagens eleitorais; \
                o p-valor não controla o erro tipo I (use null = :resampled)"
        end
    end

    pvalue = cfg.null === :benford ? p_asym :
             _resampled_pvalue(v, cfg.digit, expected, chi2, cfg, rng)

    return BenfordResult(cfg.digit, labels, counts, observed, expected,
                         n, n_excluded, chi2, df, pvalue, p_asym, cfg.null,
                         mad, _nigrini_conformity(mad, cfg.digit))
end

benford_test(x::AbstractVector{<:Integer};
             rng::AbstractRNG = Random.default_rng(), kwargs...) =
    benford_test(x, BenfordConfig(; kwargs...); rng = rng)

# Entradas não inteiras (Float64 de CSV, `missing`) passam pela coerção
# explícita de `_coerce_counts` antes de chegar ao método tipado.
benford_test(x::AbstractVector; skipmissing::Bool = false,
             rng::AbstractRNG = Random.default_rng(), kwargs...) =
    benford_test(_coerce_counts(x; skipmissing = skipmissing, name = "x"),
                 BenfordConfig(; kwargs...); rng = rng)

function Base.show(io::IO, ::MIME"text/plain", r::BenfordResult)
    printstyled(io, "Teste de Benford — $(r.digit)º dígito ($(r.digit)BL)\n";
                bold = true)
    println(io, "n = $(r.n) (excluídos: $(r.n_excluded))   null: $(r.null)")
    pmax = max(maximum(r.observed), maximum(r.expected))
    for (i, d) in enumerate(r.labels)
        print(io, "  ", d, " │", _bar(r.observed[i], pmax), " ")
        print(io, _fmt_pct(r.observed[i]))
        println(io, "  (esp. ", _fmt_pct(r.expected[i]), ")")
    end
    @printf(io, "χ²(%d) = %.2f   p = %s   MAD = %.4f\n",
            r.df, r.chi2, _fmt_p(r.pvalue), r.mad)
    if r.null === :resampled
        @printf(io, "(p assintótico contra Benford puro = %s — descritivo)\n",
                _fmt_p(r.pvalue_asymptotic))
    end
    print(io, "conformidade (Nigrini): ")
    printstyled(io, String(r.conformity);
                color = r.conformity in (:conforme, :aceitavel) ? :green :
                        r.conformity === :marginal ? :yellow : :red,
                bold = true)
    println(io)
    # MAD ignora n; χ² não. Em amostras grandes um desvio pequeno é
    # significativo, e as duas leituras divergem — isso não é inconsistência.
    if (r.pvalue < 0.05) != (r.conformity === :nao_conforme)
        printstyled(io, "nota: MAD e χ² discordam — o MAD mede o tamanho do \
                        desvio, o χ² mede a evidência contra H₀ (cresce com n)\n";
                    color = :cyan)
    end
end
