# ElectionForensics.jl

Forense eleitoral em Julia: teste de frações coarse de **Rozenas (2017)**,
lei de **Benford** (1º e 2º dígitos) e testes de **último/penúltimo dígito**
(Beber & Scacco 2012), com relatório integrado no terminal.

## Instalação

```julia
] add https://github.com/dantebertuzzi/ElectionForensics.jl
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
benford_test(votes; digit = 2)                    # 2BL, null reamostrado
last_digit_test(votes)                            # Beber & Scacco
last_digit_test(votes; position = :penultimate)   # ⚠ null uniforme frágil
rozenas_test(votes, totals; boundary = :logit, B = 999)
coarse_fractions(10)                              # 1/2, 1/3, 2/3, 1/4, ...
```

Contagens vindas de CSV (`Float64`, `missing`) são coagidas de forma explícita:
valores fracionários e `NaN` são erro, e `missing` exige `skipmissing = true`.

Os três resultados implementam Tables.jl — uma linha por dígito ou por fração:

```julia
using Tables, DataFrames
DataFrame(benford_test(votes))     # digit, count, observed, expected, excess, ...
DataFrame(rozenas_test(votes, totals))  # fraction, observed, null_mean, zscore, qvalue
```

## Métodos

**Rozenas (2017).** Conta seções cujo percentual `votes/totals` é
*exatamente* uma fração irredutível k/d com d ≤ `max_denom` (comparação em
aritmética racional, sem erro de float). A distribuição nula vem de bootstrap
paramétrico: percentuais perturbados por kernel gaussiano e contagens
reamostradas de `Binomial(totals[i], p̃ᵢ)`. P-valor unilateral com correção
`(1 + #{T* ≥ T}) / (B + 1)`.

O jitter é aplicado na escala **logit** (`boundary = :logit`, default): perto
de 0 e 1 o passo em escala de probabilidade encolhe sozinho, então a nula não
empurra a massa das seções quase unânimes para o miolo, onde as frações coarse
são densas. Com jitter em escala de probabilidade (`:reflect`), sob shares em U
o erro tipo I caía a 0,000 e o poder a ε = 2 % era 0,48 contra 0,67.

Esta é uma **aproximação frequentista** do método do artigo, que usa um
estimador bayesiano de densidade kernel reamostrada. Os excessos por fração são
exploratórios: use os `qvalues` (Benjamini–Hochberg), não os z brutos.

**Benford.** χ² de aderência + MAD com limiares de conformidade de
Nigrini (2012). Default é o **2º dígito** (2BL, Mebane 2008) com
`null = :resampled`. Contagens eleitorais limpas **não** seguem a lei de
Benford quando as seções têm tamanho homogêneo: contra o null clássico o teste
rejeita eleições simuladas limpas em 72 % dos casos (1BL: 100 %). O null
reamostrado — jitter gaussiano em `log10(contagem)` — devolve a taxa a 0,05
e responde à pergunta útil: "o dígito é anômalo *dado* o formato empírico das
contagens?". Use `null = :benford` apenas como estatística descritiva.

**Último dígito.** χ² contra uniforme em 0:9, excluindo contagens
pequenas (`min_value = 10`). Reporta também `freq{0,5}` (esperado ≈ 0.20;
números fabricados superusam 0 e 5). Calibrado para seções com ≥ 60 eleitores.
O **penúltimo** dígito exige suavidade da densidade numa escala de 100 e por
isso é opt-in: com seções de 150–900 eleitores seu null uniforme rejeita
eleições limpas em 50 % dos casos (100 % com seções de 100–200). O campo
`null_valid` sinaliza quando a condição não é atendida.

## Meça a calibração antes de concluir

A validade do null de cada teste depende do formato da sua distribuição de
contagens, e nenhuma estatística simples prediz quando ele falha. Em vez de
prometer, o pacote mede:

```julia
calibration_check(votes, totals)   # simula eleições limpas com os SEUS totais
```

Ele reporta a taxa de erro tipo I que cada teste de fato entrega nesses dados.
`anticonservador` significa que o teste rejeita eleições limpas acima de α —
não conclua fraude a partir dele. `forensics_report` roda essa checagem por
default e a imprime antes dos p-valores.

## Advertências

Nenhum teste isolado prova fraude. Falsos positivos surgem de seções
pequenas, arredondamento administrativo e agregação; leia os testes em
conjunto e no contexto institucional. O jitter em torno dos percentuais
observados torna o teste de Rozenas levemente conservador sob fraude
maciça — rejeições são evidência forte.

**O que estes testes não detectam.** Nenhum deles tem poder contra *ballot
stuffing* proporcional: em simulação, adicionar 15 % de votos a 20 % das seções
deixa os quatro testes em α. Rozenas detecta metas percentuais redondas; o
último dígito detecta arredondamento decimal. Fraude sem assinatura em dígitos
ou em frações exatas passa despercebida.

Os excessos por fração exibidos pelo teste de Rozenas são **exploratórios**:
sob H₀ o maior z entre as 31 frações tem mediana 2,4 e passa de 2 em 68 % das
eleições limpas. Use os `qvalues` (Benjamini–Hochberg), não os z brutos.

## Validação

`validation/` contém os scripts de calibração e poder que sustentam as
afirmações acima. Rode-os com `julia -t auto validation/<script>.jl`.

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
