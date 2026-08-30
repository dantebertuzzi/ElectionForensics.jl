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
- `calibration_check` — auto-diagnóstico que mede, por simulação com os
  `totals` do usuário, a taxa de erro tipo I que cada teste de fato entrega
  nesses dados, classificando cada um como `:calibrado`, `:conservador` ou
  `:anticonservador`. `forensics_report` o executa por default e imprime o
  resultado antes de qualquer p-valor.
- `boundary = :logit | :reflect` em `rozenas_test`
- `null = :kernel | :betabinomial` em `rozenas_test`, com
  `fit_betabinomial_mixture` — implementação do modelo do pacote `spikes` de
  Rozenas: mistura de Beta-Binomiais por EM e reamostragem da posterior. Duas
  divergências deliberadas em relação ao original: seleção de componentes por
  BIC em vez da correlação entre densidades de L consecutivos, e teto de
  iterações em todos os laços. Toda mistura carrega o diagnóstico de aderência
  `misfit`
- Interface Tables.jl nos cinco tipos de resultado
- Coerção explícita de entrada: `Float64` com valor inteiro é aceito; valores
  fracionários e `NaN` são erro; `missing` exige `skipmissing = true`
- `election_fingerprint` e `log_vote_rate` — fingerprint de Klimek et al.
  (2012): histograma 2-D comparecimento × votos no vencedor, e assimetria e
  curtose da taxa logarítmica de voto. Diagnóstico visual, sem p-valor. O
  modelo paramétrico de fraude do artigo não está implementado (especificação
  no Supporting Information, fora do preprint)
- `FingerprintConfig` — `bins` (resolução do histograma) e
  `vote_axis = :electorate | :valid`, que escolhe o denominador do eixo
  vertical. O default `:electorate` (`Wᵢ/Nᵢ`) é o do artigo; `:valid`
  (`Wᵢ/Vᵢ`) é a convenção de boa parte das reimplementações e muda a leitura
  do eixo
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

### Validação em dados reais
- 29 pleitos americanos de 2020 por precinct (OpenElections, transcritos de
  boletins oficiais estaduais): o null clássico de Benford rejeita **8/29**
  eleições legítimas (KS contra U(0,1): p = 0,0001); o null reamostrado
  rejeita **0/29** (KS = 0,16).
- Rússia 2012, 95 413 seções (Shpilkin; dataset de Kobak, Shpilkin &
  Pshenichnikov 2016): o teste de frações coarse detecta a anomalia
  (p = 0,043), com 4/5 e 3/4 no topo — os percentuais redondos que o artigo
  documenta.
- χ², p-valor e MAD conferidos contra `scipy` (Δ ≤ 2,8e-12) e `benford_py`
  (Δ ≤ 1e-14) nos mesmos dados reais.

### Limitações conhecidas
- O fingerprint exige o **eleitorado** por seção, que a maioria dos
  repositórios de resultados não publica — só 6 dos 40 estados do
  OpenElections o fazem.
- A nula `:betabinomial` é paramétrica: um desajuste numa faixa fina da
  distribuição de percentuais fabrica excesso nas frações que caem ali. Nos
  dados russos de 2012 ela subestima a cauda esquerda em 3,8× e produz dez
  frações com q < 0,05 entre 1/7 e 1/4. Por isso `:kernel` é o default e o
  campo `misfit` avisa acima de 2.
- Nenhum teste da bateria detecta *ballot stuffing* proporcional.
- O teste do último dígito é anti-conservador quando muitas contagens ficam
  entre 10 e 50 — chega a 0,37 com seções de 10 a 120 eleitores. Não há
  `min_value` que resolva em todos os regimes (testado em
  `validation/18_criterio_ultimo_digito.jl`), por isso a resposta do pacote é
  `calibration_check`, que mede e avisa.
- O teste de frações coarse exige percentual *exatamente* igual a k/d, o que
  requer `totals` divisível por d. Nos dados russos, 77 % das seções são
  elegíveis para alguma fração com d ≤ 10 mas só 2,6 % acertam uma — o teste vê
  apenas a fatia da manipulação que cai numa fração exata. `max_denom` maior
  recupera parte do sinal (p = 0,042 → 0,017 com d ≤ 100).
- Nenhum pleito brasileiro foi testado: o CDN do TSE bloqueia clientes
  não-browser. `validation/14_dados_reais_tse.jl` está pronto para rodar.
