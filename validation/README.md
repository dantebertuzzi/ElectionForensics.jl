# Scripts de validação

Gerados na revisão de 2026-08-29. Rode com o ambiente de auditoria:

```bash
julia --project=. -e 'using Pkg; Pkg.add(["HypothesisTests","Aqua","JET"])'
julia --project=. -t auto validation/01_calibracao_benford.jl
```

| script | pergunta | achado |
|---|---|---|
| `01_calibracao_benford.jl` | p-valores de Benford são U(0,1) sob H0? | 1BL rejeita 100 % das eleições limpas; 2BL até 25 % |
| `01b_benford_isolado.jl` | o erro está na maquinaria χ² ou no DGP? | **maquinaria correta** (α exato); o problema é o modelo |
| `02_calibracao_ultimodigito.jl` | último dígito é calibrado? | sim para seções ≥ 60 eleitores |
| `03_calibracao_rozenas.jl` | erro tipo I do teste de Rozenas | ≤ α em todos os 12 cenários |
| `04_rozenas_diagnostico.jl` | a nula reproduz a variância verdadeira? | sim (razão 0,87–0,99); `h` é robusto |
| `05_rozenas_clamp.jl` | o `clamp(·,0,1)` distorce a nula? | sim: até 22 % dos sorteios colados na fronteira |
| `06_poder.jl` | curvas de poder por tipo de fraude | penúltimo dígito: 46–52 % de falsos positivos |
| `07_penultimo_digito.jl` | caracteriza o falso positivo do penúltimo | 50 % (seções 150–900), 100 % (100–200) |
| `08_criterio_diagnostico.jl` | um guarda por amplitude salva os testes? | não para o penúltimo |
| `09_casos_limite.jl` | n=1, zeros, empates, `missing`, tipos | vários casos degenerados sem guarda |
| `10a/10b` | referência cruzada com `scipy` | concordância a 2,8e-12 |
| `11_zscores_selecao.jl` | os z por fração são interpretáveis? | não: maior z tem mediana 2,4 sob H0 |
| `12_guarda_aplicabilidade.jl` | amplitude prediz a validade da 2BL? | não de forma confiável |
| `13_benford_null_reamostrado.jl` | null reamostrado corrige a 2BL? | **sim**: 0,72 → 0,046 |
| `14_dados_reais_tse.jl` | dados reais do TSE | **não executado** (CDN do TSE: HTTP 403) |
| `fetch_dados_reais.py` | baixa EUA 2020 (OpenElections) e Rússia 2012 (Shpilkin) | — |
| `16_dados_reais.jl` | **dados reais**: calibração em 29 pleitos americanos, poder na Rússia 2012 | 2BL clássica rejeita 8/29 eleições legítimas; reamostrada 0/29. Rússia detectada (p=0,043) |
| `17_heterogeneidade_secoes.jl` | a nula do Rozenas aguenta seções de N=2 a N=5000? | sim: erro tipo I 0,011–0,055 com totais reais |
|_(o patch da primeira rodada foi aplicado; ver `git log`)_|||
| `15_rozenas_escala.jl` | qual escala de jitter calibra a nula do Rozenas? | **logit**: erro tipo I 0,000 → 0,030 e poder 0,48 → 0,67 sob shares em U |
| `18_criterio_ultimo_digito.jl` | existe diagnóstico escalar que prediga a falha do null uniforme? | **não**: correlações de +0,25 a −0,55 com contraexemplos nos dois sentidos |
| `19_valida_calibration_check.jl` | `calibration_check` acerta a taxa que prevê? | sim: 0,372 medido vs 0,364 previsto no pior regime |
| `20_concorrentes.py` | as bibliotecas existentes acusam eleições legítimas? | sim: 66–97 % dos 29 pleitos, contra 0 % do null reamostrado |
| `21_nula_betabinomial.jl` | a nula Beta-Binomial do `spikes` é melhor? | em simulação, equivalente; em dados reais, desajusta a cauda e inventa 10 frações |
| `22_fingerprint.jl` | reproduz a Fig. 3 de Klimek et al. (2012) | Illinois em (−0,01; 3,00); Rússia em (−2,15; 9,12) com 2,81 % das seções no canto |
