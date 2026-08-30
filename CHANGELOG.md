# Changelog

Formato baseado em [Keep a Changelog](https://keepachangelog.com/pt-BR/1.1.0/);
versionamento segue [SemVer](https://semver.org/lang/pt-BR/).

## [0.1.0] — não lançado

Primeira versão pública. Testes forenses de dígitos e de frações coarse, com a
calibração de cada null verificada por simulação (`validation/`).

### Adicionado
- `benford_test` (1BL/2BL) com null reamostrado, configurável por `BenfordConfig`
- `last_digit_test` (último e penúltimo dígito), configurável por `DigitTestConfig`
- `rozenas_test` — frações coarse com nula por bootstrap paramétrico
- `forensics_report` — bateria integrada com saída colorida
- `qvalues` (Benjamini–Hochberg) por fração em `RozenasResult`
- `null_valid` em `LastDigitResult`
- `boundary = :logit | :reflect` em `rozenas_test`
- Interface Tables.jl nos três tipos de resultado
- Coerção explícita de entrada: `Float64` com valor inteiro é aceito; valores
  fracionários e `NaN` são erro; `missing` exige `skipmissing = true`
- Suíte `validation/` com calibração sob H₀, curvas de poder e referência
  cruzada com `scipy`

### Notas de calibração
Três nulls clássicos da literatura **não** são válidos para contagens
eleitorais, e o pacote se afasta deles deliberadamente:

- **Lei de Benford.** Com seções de tamanho homogêneo, o χ² contra Benford
  rejeita eleições limpas em 72 % dos casos (1BL: 100 %). O default
  `null = :resampled` compara o dígito contra a distribuição empírica de
  contagens suavizada em `log10`, devolvendo a taxa a 0,054.
  `null = :benford` continua disponível como estatística descritiva.
- **Penúltimo dígito.** O null uniforme rejeita eleições limpas em até 100 %
  dos casos quando as contagens cobrem menos de duas décadas. É opt-in em
  `forensics_report` e `null_valid` sinaliza a condição.
- **Jitter do teste de Rozenas.** Na escala de probabilidade, a banda global
  empurra seções quase unânimes para o miolo: erro tipo I de 0,000 e poder de
  0,48 a ε = 2 %. Na escala logit, 0,030 e 0,70.

### Limitações conhecidas
- Nenhum teste da bateria detecta *ballot stuffing* proporcional.
- O teste do último dígito continua levemente anti-conservador (0,08 contra
  α = 0,05) quando a distribuição de percentuais é fortemente bimodal e muitas
  contagens ficam entre 10 e 30; um aviso é emitido.
- A calibração foi verificada em simulação, não contra um pleito real.
