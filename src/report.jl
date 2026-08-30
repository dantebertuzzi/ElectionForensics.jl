# ── relatório integrado ──────────────────────────────────────────────

"""
    forensics_report(votes::AbstractVector{<:Integer},
                     totals::AbstractVector{<:Integer};
                     benford_digit::Int = 2,
                     benford_null::Symbol = :resampled,
                     penultimate::Bool = false,
                     max_denom::Int = 10,
                     B::Int = 999,
                     rng::AbstractRNG = Random.default_rng(),
                     io::IO = stdout)
        -> NamedTuple

Roda a bateria de testes forenses sobre as contagens de votos de um candidato
(`votes`) e os totais por seção (`totals`), imprimindo os resultados no terminal:

1. Benford `benford_digit`-ésimo dígito sobre `votes` (null reamostrado)
2. Último dígito de `votes` (Beber & Scacco)
3. Frações coarse de Rozenas (2017) sobre `votes ./ totals`
4. *(opcional)* Penúltimo dígito, se `penultimate = true`

Retorna `(benford = ..., last_digit = ..., penultimate = ..., rozenas = ...)`;
`penultimate` é `nothing` quando não solicitado.

!!! warning "O que a bateria NÃO detecta"
    Nenhum destes testes tem poder contra *ballot stuffing* proporcional
    (adicionar uma fração constante dos votos a seções escolhidas): sob fraude
    de 20 % das seções o poder dos quatro testes fica em α
    (`validation/06_poder.jl`). O teste de Rozenas detecta metas percentuais
    redondas; o de último dígito detecta arredondamento decimal. Fraude que não
    deixa assinatura em dígitos ou em frações exatas passa despercebida.

    `penultimate = true` é desativado por default porque seu null uniforme
    rejeita eleições limpas em até 100 % dos casos — veja [`last_digit_test`](@ref).

Nenhum teste isolado prova fraude; leia os resultados em conjunto e à luz do
contexto (tamanho das seções, regras eleitorais, agregação).
"""
function forensics_report(votes::AbstractVector{<:Integer},
                          totals::AbstractVector{<:Integer};
                          benford_digit::Int = 2,
                          benford_null::Symbol = :resampled,
                          penultimate::Bool = false,
                          max_denom::Int = 10,
                          B::Int = 999,
                          rng::AbstractRNG = Random.default_rng(),
                          io::IO = stdout)
    length(votes) == length(totals) ||
        throw(DimensionMismatch("votes e totals devem ter o mesmo comprimento"))
    sum(totals) > 0 || throw(ArgumentError("sum(totals) deve ser > 0"))
    hr = "─"^60

    printstyled(io, hr, "\n"; color = :cyan)
    printstyled(io, " ElectionForensics — relatório\n"; bold = true, color = :cyan)
    printstyled(io, hr, "\n"; color = :cyan)
    println(io, "seções: ", length(votes),
            "   votos: ", sum(votes), " / ", sum(totals),
            @sprintf(" (%.1f%%)", 100 * sum(votes) / sum(totals)))
    println(io)

    # Cada teste tem um mínimo estrutural de contagem. Numa eleição municipal
    # pequena, exigi-los todos abortaria o relatório inteiro; preferimos pular
    # o teste inaplicável e dizer por quê.
    ben = _try_test(io, "Benford $(benford_digit)BL",
                    () -> benford_test(votes, BenfordConfig(digit = benford_digit,
                                       null = benford_null, B = B); rng = rng))
    ld = _try_test(io, "último dígito",
                   () -> last_digit_test(votes; position = :last))
    pd = penultimate ?
         _try_test(io, "penúltimo dígito",
                   () -> last_digit_test(votes; position = :penultimate)) : nothing

    rz = _try_test(io, "frações coarse (Rozenas)",
                   () -> rozenas_test(votes, totals;
                                      fractions = coarse_fractions(max_denom),
                                      B = B, rng = rng))
    printstyled(io, hr, "\n"; color = :cyan)

    return (benford = ben, last_digit = ld, penultimate = pd, rozenas = rz)
end

# Roda um teste da bateria, imprime o resultado e devolve-o; se o teste for
# inaplicável aos dados (contagens abaixo do mínimo estrutural), registra o
# motivo e devolve `nothing` em vez de abortar o relatório.
function _try_test(io::IO, nome::AbstractString, f)
    try
        r = f()
        show(io, MIME"text/plain"(), r)
        println(io)
        return r
    catch e
        e isa ArgumentError || rethrow()
        printstyled(io, "· $nome: não aplicável — ", e.msg, "\n\n"; color = :yellow)
        return nothing
    end
end
