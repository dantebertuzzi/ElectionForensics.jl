"""
    ElectionForensics

Forense eleitoral em Julia: testes de dígitos (Benford 1BL/2BL, último
dígito de Beber & Scacco) e o teste de frações coarse de Rozenas (2017)
para detectar manipulação de percentuais de votos.

Funções principais:

- [`benford_test`](@ref) — conformidade com a lei de Benford (1º ou 2º dígito)
- [`last_digit_test`](@ref) — uniformidade do último/penúltimo dígito
- [`rozenas_test`](@ref) — excesso de precintos em frações "redondas" (k/d)
- [`forensics_report`](@ref) — roda a bateria completa e imprime no terminal

Referências:
- Rozenas, A. (2017). "Detecting Election Fraud from Irregularities in
  Vote-Share Distributions". *Political Analysis* 25(1).
- Mebane, W. (2008). "Election Forensics: The Second-Digit Benford's Law Test".
- Beber, B. & Scacco, A. (2012). "What the Numbers Say: A Digit-Based Test
  for Election Fraud". *Political Analysis* 20(2).
- Nigrini, M. (2012). *Benford's Law*. (limiares de MAD)
- Deckert, Myagkov & Ordeshook (2011) — advertência sobre o 1º dígito
  em dados eleitorais.
"""
module ElectionForensics

using Distributions
using Statistics
using Random
using Printf

export first_digit, second_digit, last_digit, penultimate_digit,
       benford_test, last_digit_test, rozenas_test, coarse_fractions,
       forensics_report,
       BenfordResult, LastDigitResult, RozenasResult

include("utils.jl")
include("benford.jl")
include("lastdigit.jl")
include("rozenas.jl")
include("report.jl")

end # module
