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
```

## Configuração

```@docs
BenfordConfig
DigitTestConfig
CalibrationConfig
```

## Resultados

Os três tipos de resultado implementam a interface de Tables.jl: uma linha por
dígito (Benford, último dígito) ou por fração (Rozenas).

```@docs
BenfordResult
LastDigitResult
RozenasResult
CalibrationResult
```

## Utilidades

```@docs
coarse_fractions
first_digit
second_digit
last_digit
penultimate_digit
```
