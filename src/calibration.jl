# ── auto-diagnóstico de calibração ───────────────────────────────────
#
# Os nulls destes testes valem sob condições que dependem do FORMATO dos dados,
# e nenhuma estatística escalar simples prediz quando falham (testamos várias em
# validation/18_criterio_ultimo_digito.jl: as correlações vão de +0,25 a −0,55,
# com contraexemplos em ambas as direções).
#
# Em vez de prometer uma garantia que não existe, o pacote mede: simula
# eleições limpas com os SEUS tamanhos de seção e um formato de percentuais
# como o seu, roda o teste, e reporta a taxa de rejeição que ele de fato
# entrega. Se sair muito acima de α, o p-valor do teste não é confiável
# NESSES dados — e o usuário fica sabendo antes de concluir qualquer coisa.

"""
    CalibrationConfig

Parâmetros do [`calibration_check`](@ref).

- `alpha`: nível nominal a verificar.
- `R`: eleições limpas simuladas (precisão da taxa ≈ `sqrt(α(1-α)/R)`).
- `B`: réplicas internas dos testes de Monte Carlo.
- `bandwidth`: banda do kernel usado para suavizar os percentuais observados
  antes de reamostrar. `nothing` usa Silverman na escala logit.
- `tests`: quais testes verificar.
"""
Base.@kwdef struct CalibrationConfig
    alpha::Float64 = 0.05
    R::Int = 200
    B::Int = 199
    bandwidth::Union{Nothing,Float64} = nothing
    tests::Vector{Symbol} = [:benford, :last_digit, :penultimate, :rozenas]
end

"""
    CalibrationResult

Resultado de [`calibration_check`](@ref). Implementa a interface de Tables.jl:
uma linha por teste, com `test`, `alpha`, `rejection_rate`, `se` (erro padrão
de Monte Carlo), `ratio` (taxa ÷ α) e `verdict`.

`verdict` é `:calibrado` quando a taxa observada é compatível com α dentro de
dois erros padrão, `:anticonservador` quando fica significativamente acima
(o teste acusa fraude em dados limpos) e `:conservador` quando fica abaixo
(o teste perde poder).
"""
struct CalibrationResult
    tests::Vector{Symbol}
    alpha::Float64
    rejection_rate::Vector{Float64}
    se::Vector{Float64}
    verdict::Vector{Symbol}
    R::Int
    m::Int
end

function _verdict(rate::Float64, alpha::Float64, se::Float64)
    se ≤ 0 && return :calibrado
    rate > alpha + 2se && return :anticonservador
    rate < alpha - 2se && return :conservador
    return :calibrado
end

"""
    calibration_check(votes, totals; kwargs...) -> CalibrationResult
    calibration_check(votes, totals, cfg::CalibrationConfig; rng) -> CalibrationResult

Mede a taxa de erro tipo I que cada teste **de fato** entrega em dados com o
formato dos seus.

Simula `R` eleições limpas por construção: os percentuais latentes vêm de uma
reamostragem suavizada dos percentuais observados (kernel gaussiano na escala
logit, que preserva o formato sem herdar a estrutura de dígitos), e as
contagens de `Binomial(totals[i], p̃ᵢ)` com os **seus** `totals`. Como H₀ é
verdadeiro nessas réplicas, a fração de rejeições a `alpha` estima o erro
tipo I real.

Interpretação:

- `:calibrado` — o p-valor do teste significa o que promete nestes dados.
- `:anticonservador` — o teste rejeita eleições limpas com frequência acima de
  `alpha`. **Não conclua fraude a partir dele nestes dados.**
- `:conservador` — o teste perde poder; não rejeitar não é evidência de nada.

Rode isto antes de interpretar [`forensics_report`](@ref). O custo é
`R × B` reamostragens: com os defaults, alguns segundos para milhares de seções.

!!! note "Por que isso existe"
    A validade do null uniforme (último e penúltimo dígito) e do null de
    Benford depende do formato da distribuição de contagens, e não há
    estatística escalar simples que prediga a falha — várias foram testadas em
    `validation/18_criterio_ultimo_digito.jl`. Medir é mais honesto que
    prometer.

# Exemplo

```julia
cal = calibration_check(votes, totals)
using Tables; Tables.columntable(cal)
```
"""
function calibration_check(votes::AbstractVector{<:Integer},
                           totals::AbstractVector{<:Integer},
                           cfg::CalibrationConfig;
                           rng::AbstractRNG = Random.default_rng())
    m = length(votes)
    m == length(totals) ||
        throw(DimensionMismatch("votes e totals devem ter o mesmo comprimento"))
    m > 0 || throw(ArgumentError("vetores vazios"))
    0 < cfg.alpha < 1 || throw(ArgumentError("alpha deve estar em (0,1)"))
    cfg.R ≥ 20 || throw(ArgumentError("R deve ser ≥ 20"))
    isempty(cfg.tests) && throw(ArgumentError("nenhum teste selecionado"))
    for t in cfg.tests
        t in (:benford, :last_digit, :penultimate, :rozenas) ||
            throw(ArgumentError("teste desconhecido: $t"))
    end
    @inbounds for i in 1:m
        totals[i] > 0 || throw(ArgumentError("totals[$i] deve ser > 0"))
        0 ≤ votes[i] ≤ totals[i] ||
            throw(ArgumentError("votes[$i] fora de [0, totals[$i]]"))
    end

    T = [Int(x) for x in totals]
    jshares = [_compress(Float64(votes[i]), Float64(T[i])) for i in 1:m]
    lz = _logit.(jshares)
    h = cfg.bandwidth === nothing ? min(first(_silverman(lz)), _HLOGIT_CAP) :
                                    Float64(cfg.bandwidth)

    contagem = Dict(t => 0 for t in cfg.tests)
    validas  = Dict(t => 0 for t in cfg.tests)
    v = Vector{Int}(undef, m)
    for _ in 1:cfg.R
        @inbounds for i in 1:m
            p = _ilogit(lz[rand(rng, 1:m)] + h * randn(rng))
            v[i] = rand(rng, Binomial(T[i], p))
        end
        for t in cfg.tests
            p = try
                if t === :benford
                    benford_test(v, BenfordConfig(B = cfg.B, warn = false); rng = rng).pvalue
                elseif t === :last_digit
                    last_digit_test(v, DigitTestConfig(position = :last, warn = false)).pvalue
                elseif t === :penultimate
                    last_digit_test(v, DigitTestConfig(position = :penultimate, warn = false)).pvalue
                else
                    rozenas_test(v, T; B = cfg.B, rng = rng).pvalue
                end
            catch e
                e isa ArgumentError || rethrow()
                continue          # teste inaplicável nesta réplica
            end
            validas[t] += 1
            p < cfg.alpha && (contagem[t] += 1)
        end
    end

    taxas = Float64[]; ses = Float64[]; vers = Symbol[]
    for t in cfg.tests
        n = validas[t]
        taxa = n > 0 ? contagem[t] / n : NaN
        se = n > 0 ? sqrt(cfg.alpha * (1 - cfg.alpha) / n) : NaN
        push!(taxas, taxa); push!(ses, se)
        push!(vers, isnan(taxa) ? :inaplicavel : _verdict(taxa, cfg.alpha, se))
    end
    return CalibrationResult(copy(cfg.tests), cfg.alpha, taxas, ses, vers, cfg.R, m)
end

calibration_check(votes::AbstractVector{<:Integer}, totals::AbstractVector{<:Integer};
                  rng::AbstractRNG = Random.default_rng(), kwargs...) =
    calibration_check(votes, totals, CalibrationConfig(; kwargs...); rng = rng)

function Base.show(io::IO, ::MIME"text/plain", r::CalibrationResult)
    printstyled(io, "Auto-diagnóstico de calibração\n"; bold = true)
    @printf(io, "%d seções · %d eleições limpas simuladas · α = %.3f\n",
            r.m, r.R, r.alpha)
    println(io, "  teste            rejeição   esperado    veredito")
    nomes = Dict(:benford => "Benford", :last_digit => "último dígito",
                 :penultimate => "penúltimo dígito", :rozenas => "frações coarse")
    for (i, t) in enumerate(r.tests)
        @printf(io, "  %-16s %-10s %-11s ", nomes[t],
                isnan(r.rejection_rate[i]) ? "—" : @sprintf("%.3f", r.rejection_rate[i]),
                @sprintf("%.3f±%.3f", r.alpha, 2r.se[i]))
        printstyled(io, String(r.verdict[i]);
                    color = r.verdict[i] === :anticonservador ? :red :
                            r.verdict[i] === :conservador ? :yellow : :green,
                    bold = r.verdict[i] === :anticonservador)
        println(io)
    end
    if any(==(:anticonservador), r.verdict)
        printstyled(io, "⚠ testes anticonservadores rejeitam dados LIMPOS acima de α; \
                        não conclua fraude a partir deles nestes dados\n"; color = :red)
    end
end
