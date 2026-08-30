# Referência

```@meta
CurrentModule = ElectionForensics
```

## Testes

```@docs
benford_test
last_digit_test
rozenas_test
forensics_report
calibration_check
fit_betabinomial_mixture
election_fingerprint
log_vote_rate
```

## Configuração

```@docs
BenfordConfig
DigitTestConfig
CalibrationConfig
FingerprintConfig
```

## Resultados

Os três tipos de resultado implementam a interface de Tables.jl: uma linha por
dígito (Benford, último dígito) ou por fração (Rozenas).

```@docs
BenfordResult
LastDigitResult
RozenasResult
CalibrationResult
BetaBinomialMixture
FingerprintResult
```

## Utilidades

```@docs
coarse_fractions
first_digit
second_digit
last_digit
penultimate_digit
```
