# ── relatório integrado ──────────────────────────────────────────────

"""
    forensics_report(votes::AbstractVector{<:Integer},
                     totals::AbstractVector{<:Integer};
                     benford_digit::Int = 2,
                     max_denom::Int = 10,
                     B::Int = 999,
                     rng::AbstractRNG = Random.default_rng(),
                     io::IO = stdout)
        -> NamedTuple

Roda a bateria completa de testes forenses sobre as contagens de votos
de um candidato (`votes`) e os totais por seção (`totals`), imprimindo
os resultados no terminal:

1. Benford `benford_digit`-ésimo dígito sobre `votes`
2. Último dígito de `votes` (Beber & Scacco)
3. Penúltimo dígito de `votes`
4. Frações coarse de Rozenas (2017) sobre `votes ./ totals`

Retorna `(benford = ..., last_digit = ..., penultimate = ..., rozenas = ...)`
para inspeção posterior.

Nenhum teste isolado prova fraude; leia os quatro em conjunto e à luz do
contexto (tamanho das seções, regras eleitorais, agregação).
"""
function forensics_report(votes::AbstractVector{<:Integer},
                          totals::AbstractVector{<:Integer};
                          benford_digit::Int = 2,
                          max_denom::Int = 10,
                          B::Int = 999,
                          rng::AbstractRNG = Random.default_rng(),
                          io::IO = stdout)
    hr = "─"^60

    printstyled(io, hr, "\n"; color = :cyan)
    printstyled(io, " ElectionForensics — relatório\n"; bold = true, color = :cyan)
    printstyled(io, hr, "\n"; color = :cyan)
    println(io, "seções: ", length(votes),
            "   votos: ", sum(votes), " / ", sum(totals),
            @sprintf(" (%.1f%%)", 100 * sum(votes) / sum(totals)))
    println(io)

    ben = benford_test(votes; digit = benford_digit)
    show(io, MIME"text/plain"(), ben)
    println(io)

    ld = last_digit_test(votes; position = :last)
    show(io, MIME"text/plain"(), ld)
    println(io)

    pd = last_digit_test(votes; position = :penultimate)
    show(io, MIME"text/plain"(), pd)
    println(io)

    rz = rozenas_test(votes, totals;
                      fractions = coarse_fractions(max_denom), B = B, rng = rng)
    show(io, MIME"text/plain"(), rz)
    printstyled(io, hr, "\n"; color = :cyan)

    return (benford = ben, last_digit = ld, penultimate = pd, rozenas = rz)
end
