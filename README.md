# ElectionForensics.jl

Forense eleitoral em Julia: teste de frações coarse de **Rozenas (2017)**,
lei de **Benford** (1º e 2º dígitos) e testes de **último/penúltimo dígito**
(Beber & Scacco 2012), com relatório integrado no terminal.

## Instalação

```julia
] dev /caminho/para/ElectionForensics.jl
```

## Uso rápido

```julia
using ElectionForensics

# votes[i]  = votos do candidato na seção i
# totals[i] = votos válidos na seção i

res = forensics_report(votes, totals)   # bateria completa, saída colorida

res.rozenas.pvalue
res.benford.conformity
```

Testes individuais:

```julia
benford_test(votes; digit = 2)                    # 2BL (Mebane)
last_digit_test(votes)                            # Beber & Scacco
last_digit_test(votes; position = :penultimate)
rozenas_test(votes, totals; max_denom = 10, B = 999)
coarse_fractions(10)                              # 1/2, 1/3, 2/3, 1/4, ...
```

## Métodos

**Rozenas (2017).** Conta seções cujo percentual `votes/totals` é
*exatamente* uma fração irredutível k/d com d ≤ `max_denom` (comparação em
`Rational`, sem erro de float). A distribuição nula vem de bootstrap
paramétrico: percentuais perturbados por kernel gaussiano (banda de
Silverman) e contagens reamostradas de `Binomial(totals[i], p̃ᵢ)` —
aproximação da densidade kernel reamostrada do artigo. P-valor unilateral
com correção `(1 + #{T* ≥ T}) / (B + 1)`.

**Benford.** χ² de aderência + MAD com limiares de conformidade de
Nigrini (2012). Default é o **2º dígito** (2BL, Mebane 2008): a aderência
do 1º dígito em contagens de votos é contestada (Deckert, Myagkov &
Ordeshook 2011).

**Último dígito.** χ² contra uniforme em 0:9, excluindo contagens
pequenas (`min_value = 10`). Reporta também `freq{0,5}` (esperado ≈ 0.20;
números fabricados superusam 0 e 5).

## Advertências

Nenhum teste isolado prova fraude. Falsos positivos surgem de seções
pequenas, arredondamento administrativo e agregação; leia os testes em
conjunto e no contexto institucional. O jitter em torno dos percentuais
observados torna o teste de Rozenas levemente conservador sob fraude
maciça — rejeições são evidência forte.

## Referências

- Rozenas, A. (2017). Detecting Election Fraud from Irregularities in
  Vote-Share Distributions. *Political Analysis* 25(1), 41–56.
- Beber, B.; Scacco, A. (2012). What the Numbers Say: A Digit-Based Test
  for Election Fraud. *Political Analysis* 20(2), 211–234.
- Mebane, W. (2008). Election Forensics: The Second-Digit Benford's Law
  Test and Recent American Presidential Elections.
- Nigrini, M. (2012). *Benford's Law*. Wiley.
- Deckert, J.; Myagkov, M.; Ordeshook, P. (2011). Benford's Law and the
  Detection of Election Fraud. *Political Analysis* 19(3), 245–268.
