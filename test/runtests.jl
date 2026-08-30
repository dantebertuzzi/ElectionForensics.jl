using ElectionForensics
using Test
using Aqua
using Tables
using Random
using Distributions
using SpecialFunctions
using Statistics

@testset "ElectionForensics.jl" begin

    @testset "extração de dígitos" begin
        @test first_digit(7) == 7
        @test first_digit(1234) == 1
        @test first_digit(987654) == 9
        @test_throws ArgumentError first_digit(0)

        @test second_digit(10) == 0
        @test second_digit(1234) == 2
        @test second_digit(95) == 5
        @test_throws ArgumentError second_digit(9)

        @test last_digit(1230) == 0
        @test last_digit(7) == 7
        @test last_digit(-7) == 7
        @test last_digit(0) == 0

        @test penultimate_digit(1234) == 3
        @test penultimate_digit(10) == 1
        @test penultimate_digit(-1234) == 3
        @test_throws ArgumentError penultimate_digit(5)
        @test_throws ArgumentError penultimate_digit(0)
    end

    @testset "coarse_fractions" begin
        fr = coarse_fractions(4)
        @test fr == sort([1//2, 1//3, 2//3, 1//4, 3//4])
        @test all(f -> 0 < f < 1, coarse_fractions(10))
        @test allunique(coarse_fractions(10))
        @test_throws ArgumentError coarse_fractions(1)
    end

    # Uma semente por testset. Com um RNG encadeado, inserir um teste no meio
    # muda os dados de todos os posteriores e falhas viram mistério.
    seed(nome) = MersenneTwister(hash(("ElectionForensics", nome)))

    @testset "Benford — dados log-uniformes conformes" begin
        rng = seed("Benford — dados log-uniformes conformes")
        # log-uniforme em [1, 10^5) segue Benford exatamente
        x = [round(Int, 10.0^(5rand(rng))) for _ in 1:5000]
        x = [xi for xi in x if xi ≥ 10]
        r1 = benford_test(x; digit = 1, null = :benford, warn = false)
        r2 = benford_test(x; digit = 2, null = :benford, warn = false)
        @test r1.pvalue > 0.001
        @test r2.pvalue > 0.001
        @test r1.conformity in (:conforme, :aceitavel)
        @test sum(r1.counts) == r1.n
        @test isapprox(sum(r1.expected), 1.0; atol = 1e-10)
        @test isapprox(sum(r2.expected), 1.0; atol = 1e-10)
    end

    @testset "Benford — dados fabricados não conformes" begin
        rng = seed("Benford — dados fabricados não conformes")
        # primeiro dígito uniforme em 1:9 viola fortemente a 1BL
        x = [d * 10^rand(rng, 1:4) + rand(rng, 0:9) for d in rand(rng, 1:9, 3000)]
        r = benford_test(x; digit = 1, null = :benford, warn = false)
        @test r.pvalue < 1e-6
        @test r.conformity == :nao_conforme
    end

    @testset "último dígito — limpo vs. fabricado" begin
        rng = seed("último dígito — limpo vs. fabricado")
        limpo = rand(rng, 100:99999, 5000)
        rl = last_digit_test(limpo)
        @test rl.pvalue > 0.001
        @test 0.15 < rl.freq_0_5 < 0.25

        # fabricado: 60% das contagens terminam em 0 ou 5
        fab = copy(limpo)
        for i in 1:3000
            fab[i] = (fab[i] ÷ 10) * 10 + rand(rng, (0, 5))
        end
        rf = last_digit_test(fab)
        @test rf.pvalue < 1e-6
        @test rf.freq_0_5 > 0.4

        rp = last_digit_test(limpo; position = :penultimate)
        @test rp.pvalue > 0.001
    end

    @testset "Rozenas — dados limpos não rejeitam" begin
        rng = seed("Rozenas — dados limpos não rejeitam")
        m = 2000
        totals = rand(rng, 150:900, m)
        p = rand(rng, Beta(8, 6), m)
        votes = [rand(rng, Binomial(totals[i], p[i])) for i in 1:m]
        r = rozenas_test(votes, totals; B = 199, rng = rng)
        @test r.pvalue > 0.05
        @test r.total_observed == sum(r.observed)
        @test length(r.null_totals) == 199
    end

    @testset "Rozenas — fraude injetada rejeita" begin
        rng = seed("Rozenas — fraude injetada rejeita")
        m = 2000
        totals = [rand(rng, 15:90) * 10 for _ in 1:m]  # divisíveis por 10
        p = rand(rng, Beta(8, 6), m)
        votes = [rand(rng, Binomial(totals[i], p[i])) for i in 1:m]
        # 8% das seções forçadas a percentuais exatos "redondos"
        alvo = (1//2, 3//5, 7//10, 3//4)
        for i in 1:round(Int, 0.08m)
            f = rand(rng, alvo)
            votes[i] = Int(totals[i] * numerator(f) ÷ denominator(f))
        end
        r = rozenas_test(votes, totals; B = 199, rng = rng)
        @test r.pvalue < 0.05
        @test r.zscore > 2
    end

    @testset "Benford — argumentos inválidos" begin
        @test_throws ArgumentError benford_test([1, 2]; digit = 3)
        @test_throws ArgumentError benford_test([10, 20]; null = :bogus)
        x = [1, 2, 3, 4, 5, 6, 7, 8, 9]  # todos < 10 para 2BL
        @test_throws ArgumentError benford_test(x; digit = 2, warn = false)
        x = [0, 0, 0]  # todos < 1 para 1BL
        @test_throws ArgumentError benford_test(x; digit = 1, warn = false)
    end

    @testset "Benford — n_excluded" begin
        x = [1, 5, 9, 10, 20, 30]  # 1BL exclui nada, 2BL exclui 1,5,9
        r1 = benford_test(x; digit = 1, warn = false)
        @test r1.n_excluded == 0
        @test r1.n == 6
        r2 = benford_test(x; digit = 2, warn = false)
        @test r2.n_excluded == 3
        @test r2.n == 3
    end

    @testset "benford_expected" begin
        e1 = ElectionForensics.benford_expected(1)
        @test length(e1) == 9
        @test isapprox(sum(e1), 1.0; atol = 1e-10)
        @test e1[1] > e1[9]  # dígito 1 mais provável que 9
        e2 = ElectionForensics.benford_expected(2)
        @test length(e2) == 10
        @test isapprox(sum(e2), 1.0; atol = 1e-10)
        @test_throws ArgumentError ElectionForensics.benford_expected(0)
        @test_throws ArgumentError ElectionForensics.benford_expected(3)
    end

    @testset "último dígito — argumentos inválidos" begin
        @test_throws ArgumentError last_digit_test([1, 2]; position = :invalid)
        x = [1, 2, 3]  # todos < 10 (default para :last)
        @test_throws ArgumentError last_digit_test(x)
        x = [10, 20, 30]
        @test_throws ArgumentError last_digit_test(x; position = :penultimate)
    end

    @testset "último dígito — min_value customizado" begin
        x = [5, 8, 9, 15, 20, 50]
        r = last_digit_test(x; position = :last, min_value = 10, warn = false)
        @test r.n == 3
        @test r.n_excluded == 3
        r2 = last_digit_test(x; position = :last, min_value = 20, warn = false)
        @test r2.n == 2
        @test r2.n_excluded == 4
        # abaixo do mínimo estrutural o null uniforme é impossível
        @test_throws ArgumentError last_digit_test(x; min_value = 5)
        @test_throws ArgumentError last_digit_test(x; position = :penultimate,
                                                   min_value = 50)
    end

    @testset "Rozenas — validação de entradas" begin
        @test_throws DimensionMismatch rozenas_test([1, 2], [10])
        @test_throws ArgumentError rozenas_test([5], [0])
        @test_throws ArgumentError rozenas_test([11], [10])
        @test_throws ArgumentError rozenas_test(Int[], Int[])
        @test_throws ArgumentError rozenas_test([10], [100]; B = 50)
    end

    @testset "Rozenas — frações customizadas e h manual" begin
        rng = seed("Rozenas — frações customizadas e h manual")
        m = 500
        totals = rand(rng, 150:900, m)
        p = rand(rng, Beta(8, 6), m)
        votes = [rand(rng, Binomial(totals[i], p[i])) for i in 1:m]
        custom_fr = [1//2, 2//3, 1//5]
        r = rozenas_test(votes, totals; fractions = custom_fr, B = 199, rng = rng)
        @test r.fractions == custom_fr
        @test r.pvalue ≥ 0
        @test length(r.null_totals) == 199

        r2 = rozenas_test(votes, totals; h = 0.05, B = 199, rng = rng)
        @test r2.h == 0.05
        @test r2.pvalue ≥ 0
    end

    @testset "Rozenas — h inválido" begin
        @test_throws ArgumentError rozenas_test([10], [100]; h = 0.0)
        @test_throws ArgumentError rozenas_test([10], [100]; h = -0.1)
    end

    @testset "forensics_report" begin
        rng = seed("forensics_report")
        m = 800
        totals = rand(rng, 100:600, m)
        votes = [rand(rng, Binomial(totals[i], 0.55)) for i in 1:m]
        res = forensics_report(votes, totals; B = 99, rng = rng,
                               io = IOBuffer())
        @test res.benford isa BenfordResult
        @test res.last_digit isa LastDigitResult
        @test res.penultimate === nothing          # opt-in desde a correção
        @test res.rozenas isa RozenasResult

        res2 = forensics_report(votes, totals; B = 99, rng = rng,
                                penultimate = true, io = IOBuffer())
        @test res2.penultimate isa LastDigitResult
    end

    @testset "Rozenas — q-valores de Benjamini–Hochberg" begin
        rng = seed("Rozenas — q-valores de Benjamini–Hochberg")
        m = 800
        totals = rand(rng, 150:900, m)
        p = rand(rng, Beta(8, 6), m)
        votes = [rand(rng, Binomial(totals[i], p[i])) for i in 1:m]
        r = rozenas_test(votes, totals; B = 199, rng = rng)
        @test length(r.qvalues) == length(r.fractions)
        @test all(0 .< r.qvalues .<= 1)
        # BH é monótono na ordem dos p brutos
        @test issorted(sort(r.qvalues))
    end

    @testset "_bh_adjust — valores analíticos" begin
        bh = ElectionForensics._bh_adjust
        @test bh(Float64[]) == Float64[]
        @test bh([0.5]) ≈ [0.5]
        # p = [.01,.02,.03,.04] , n=4 -> n*p/k = [.04,.04,.04,.04]
        @test bh([0.01, 0.02, 0.03, 0.04]) ≈ fill(0.04, 4)
        # monotonicidade forçada: p=[.01,.5] -> [.02,.5]
        @test bh([0.01, 0.5]) ≈ [0.02, 0.5]
        @test all(bh([0.9, 0.95]) .≈ 0.95)
    end

    @testset "penúltimo dígito — null_valid sinaliza inaplicabilidade" begin
        rng = seed("penúltimo dígito — null_valid sinaliza inaplicabilidade")
        # contagens numa única década: null uniforme não é defensável
        estreito = rand(rng, 150:900, 3000)
        r = last_digit_test(estreito; position = :penultimate, warn = false)
        @test r.null_valid == false
        # contagens cobrindo > 2 décadas (log-uniforme): aplicável
        largo = [round(Int, 10.0^(2 + 3rand(rng))) for _ in 1:3000]
        @test last_digit_test(largo; position = :penultimate, warn = false).null_valid
        # :last nunca é marcado inválido
        @test last_digit_test(estreito; warn = false).null_valid
    end

    @testset "Benford — null reamostrado é calibrado onde o clássico não é" begin
        # Seções homogêneas (caso TSE). A afirmação é sobre a TAXA de rejeição
        # sob H0, não sobre um sorteio: um único p-valor cruza α em ~5 % das
        # sementes por construção. 40 réplicas independentes, seeds fixas.
        nrep = 40
        rej_classico = 0
        rej_reamostr = 0
        for k in 1:nrep
            r = MersenneTwister(9000 + k)
            totals = rand(r, 200:400, 2000)
            p = rand(r, Beta(8, 6), 2000)
            votes = [rand(r, Binomial(totals[i], p[i])) for i in 1:2000]
            rej_classico += benford_test(votes; digit = 2, null = :benford,
                                         warn = false).pvalue < 0.05
            rej_reamostr += benford_test(votes; digit = 2, null = :resampled,
                                         B = 299, rng = r).pvalue < 0.05
        end
        # null clássico: falso positivo maciço (medido em ~0,72)
        @test rej_classico / nrep > 0.4
        # null reamostrado: taxa compatível com α = 0,05
        # (Binomial(40, 0.05): P(X ≥ 6) ≈ 0,016 — folga suficiente contra flake)
        @test rej_reamostr <= 5

        r0 = MersenneTwister(4242)
        totals = rand(r0, 200:400, 2000)
        p = rand(r0, Beta(8, 6), 2000)
        votes = [rand(r0, Binomial(totals[i], p[i])) for i in 1:2000]
        reamostr = benford_test(votes; digit = 2, null = :resampled, B = 299, rng = r0)
        classico = benford_test(votes; digit = 2, null = :benford, warn = false)
        @test reamostr.pvalue_asymptotic == classico.pvalue
        @test reamostr.null === :resampled
        @test reamostr.chi2 == classico.chi2
    end

    @testset "Aqua" begin
        Aqua.test_all(ElectionForensics)
    end

    @testset "coerção de entrada — Float, missing, NaN" begin
        # Float com valor inteiro é aceito
        @test benford_test([10.0, 20.0, 30.0, 45.0]; warn = false).n == 4
        @test last_digit_test([10.0, 25.0, 33.0]; warn = false).n == 3
        @test rozenas_test([50.0, 60.0], [100.0, 120.0]; B = 99).m == 2

        # Float fracionário é erro: contagem de votos não é fracionária
        @test_throws ArgumentError benford_test([10.5, 20.0]; warn = false)
        @test_throws ArgumentError benford_test([10.0, NaN]; warn = false)
        @test_throws ArgumentError benford_test([10.0, Inf]; warn = false)

        # `missing` exige consentimento explícito
        x = Union{Missing,Int}[10, 20, 30, missing]
        @test_throws ArgumentError benford_test(x; warn = false)
        @test benford_test(x; skipmissing = true, warn = false).n == 3
        @test last_digit_test(x; skipmissing = true, warn = false).n == 3
        # em rozenas descartar desalinharia os pares
        @test_throws ArgumentError rozenas_test(x, Union{Missing,Int}[1, 2, 3, 4];
                                                skipmissing = true, B = 99)

        @test_throws ArgumentError benford_test(Any["a", 10]; warn = false)
    end

    @testset "Rozenas — escala do jitter" begin
        rng = seed("Rozenas — escala do jitter")
        m = 800
        totals = rand(rng, 150:900, m)
        p = rand(rng, Beta(0.5, 0.5), m)
        votes = [rand(rng, Binomial(totals[i], p[i])) for i in 1:m]

        rl = rozenas_test(votes, totals; B = 199, rng = MersenneTwister(1))
        rr = rozenas_test(votes, totals; B = 199, boundary = :reflect,
                          rng = MersenneTwister(1))
        @test rl.boundary === :logit          # default
        @test rr.boundary === :reflect
        @test rl.h != rr.h                    # bandas em escalas diferentes
        @test 0 < rl.pvalue ≤ 1
        @test 0 < rr.pvalue ≤ 1
        @test_throws ArgumentError rozenas_test(votes, totals; boundary = :bogus)

        # Sob shares em U, a nula com reflexão superestima T e o teste fica
        # conservador; a escala logit encosta na taxa nominal.
        # (validation/15_rozenas_escala.jl mede: 0,012 vs 0,036)
        @test mean(rr.null_totals) ≥ mean(rl.null_totals)
    end

    @testset "_silverman — banda degenerada é sinalizada" begin
        h, deg = ElectionForensics._silverman([0.5, 0.5, 0.5, 0.5])
        @test deg
        @test h == ElectionForensics._H_FLOOR
        h2, deg2 = ElectionForensics._silverman(collect(0.0:0.01:1.0))
        @test !deg2
        @test h2 > 0
        # valor analítico: 0.9 · min(σ, IQR/1.34) · m^(-1/5)
        v = collect(1.0:100.0)
        h3, _ = ElectionForensics._silverman(v)
        esperado = 0.9 * min(std(v), (quantile(v, .75) - quantile(v, .25)) / 1.34) *
                   100^(-1/5)
        @test h3 ≈ esperado
    end

    @testset "conformidade de Nigrini — fronteiras fechadas" begin
        nc = ElectionForensics._nigrini_conformity
        @test nc(0.0060, 1) === :conforme       # limite superior é inclusivo
        @test nc(0.0061, 1) === :aceitavel
        @test nc(0.0120, 1) === :aceitavel
        @test nc(0.0150, 1) === :marginal
        @test nc(0.0151, 1) === :nao_conforme
        @test nc(0.0080, 2) === :conforme
        @test nc(0.0120, 2) === :marginal
        @test nc(0.0121, 2) === :nao_conforme
    end

    @testset "Tables.jl" begin
        rng = seed("Tables.jl")
        m = 600
        totals = rand(rng, 150:900, m)
        p = rand(rng, Beta(8, 6), m)
        votes = [rand(rng, Binomial(totals[i], p[i])) for i in 1:m]

        b = benford_test(votes; B = 99, rng = rng, warn = false)
        l = last_digit_test(votes; warn = false)
        z = rozenas_test(votes, totals; B = 99, rng = rng)

        for r in (b, l, z)
            @test Tables.istable(typeof(r))
            @test Tables.rowaccess(typeof(r))
            sch = Tables.schema(r)
            rows = Tables.rows(r)
            @test length(sch.names) == length(first(rows))
            @test keys(first(rows)) == sch.names
            @test length(Tables.columntable(r)[first(sch.names)]) == length(rows)
        end

        @test length(Tables.rows(b)) == 10
        @test length(Tables.rows(l)) == 10
        @test length(Tables.rows(z)) == length(z.fractions)

        # a tabela reproduz os campos do struct
        bt = Tables.columntable(b)
        @test bt.count == b.counts
        @test bt.observed ≈ b.observed
        @test sum(bt.count) == b.n
        zt = Tables.columntable(z)
        @test sum(zt.observed) == z.total_observed
        @test zt.qvalue == z.qvalues
    end

    @testset "forensics_report — degradação graciosa" begin
        # eleição minúscula: Benford e último dígito são inaplicáveis, mas o
        # relatório não pode abortar
        votes  = [3, 5, 2, 7, 4, 6, 3, 8, 5, 4, 2, 9, 6, 3, 7, 5]
        totals = [12, 14, 11, 16, 13, 15, 12, 17, 14, 13, 11, 18, 15, 12, 16, 14]
        buf = IOBuffer()
        res = forensics_report(votes, totals; B = 99, io = buf,
                               rng = MersenneTwister(5))
        @test res.benford === nothing
        @test res.last_digit === nothing
        @test res.rozenas isa RozenasResult
        saida = String(take!(buf))
        @test occursin("não aplicável", saida)
    end

    @testset "forensics_report — validação de entrada" begin
        @test_throws DimensionMismatch forensics_report([1, 2], [10];
                                                        io = IOBuffer())
        @test_throws ArgumentError forensics_report([0], [0]; io = IOBuffer())
    end

    @testset "calibration_check — validação de entrada" begin
        @test_throws DimensionMismatch calibration_check([1, 2], [10])
        @test_throws ArgumentError calibration_check(Int[], Int[])
        @test_throws ArgumentError calibration_check([5], [0])
        @test_throws ArgumentError calibration_check([11], [10])
        @test_throws ArgumentError calibration_check([5], [10]; alpha = 0.0)
        @test_throws ArgumentError calibration_check([5], [10]; alpha = 1.0)
        @test_throws ArgumentError calibration_check([5], [10]; R = 5)
        @test_throws ArgumentError calibration_check([5], [10]; tests = Symbol[])
        @test_throws ArgumentError calibration_check([5], [10]; tests = [:bogus])
    end

    @testset "calibration_check — detecta o que sabemos estar quebrado" begin
        # Regime onde o null uniforme do último dígito falha por construção:
        # seções minúsculas ⇒ contagens entre 10 e 70. A simulação de
        # validation/18 mede erro tipo I ≈ 0,37 aqui.
        r1 = MersenneTwister(31)
        tot_ruim = rand(r1, 10:120, 2000)
        p1 = rand(r1, Beta(8, 6), 2000)
        v_ruim = [rand(r1, Binomial(tot_ruim[i], p1[i])) for i in 1:2000]
        cal_ruim = calibration_check(v_ruim, tot_ruim;
                                     R = 120, tests = [:last_digit],
                                     rng = MersenneTwister(32))
        @test cal_ruim.verdict[1] === :anticonservador
        @test cal_ruim.rejection_rate[1] > 0.15

        # Regime onde o mesmo teste é sabidamente calibrado.
        r2 = MersenneTwister(33)
        tot_bom = rand(r2, 300:3000, 2000)
        p2 = rand(r2, Beta(8, 6), 2000)
        v_bom = [rand(r2, Binomial(tot_bom[i], p2[i])) for i in 1:2000]
        cal_bom = calibration_check(v_bom, tot_bom;
                                    R = 120, tests = [:last_digit],
                                    rng = MersenneTwister(34))
        @test cal_bom.verdict[1] !== :anticonservador
        @test cal_bom.rejection_rate[1] < 0.15

        # O penúltimo dígito é anticonservador em seções de tamanho típico
        # (validation/07 mede 0,503 com totais 150–900), mas fica calibrado
        # quando as contagens cobrem duas décadas — que é o que `null_valid`
        # sinaliza. O diagnóstico tem de distinguir os dois casos.
        r3 = MersenneTwister(36)
        tot_est = rand(r3, 150:900, 2000)
        p3 = rand(r3, Beta(8, 6), 2000)
        v_est = [rand(r3, Binomial(tot_est[i], p3[i])) for i in 1:2000]
        cal_pen = calibration_check(v_est, tot_est;
                                    R = 120, tests = [:penultimate],
                                    rng = MersenneTwister(37))
        @test cal_pen.verdict[1] === :anticonservador
        @test cal_pen.rejection_rate[1] > 0.2
        @test last_digit_test(v_est; position = :penultimate, warn = false).null_valid == false

        # contagens cobrindo duas décadas: o mesmo teste passa
        # (validation/07 mede 0,053 com totais uniformes em 100–9999)
        r4 = MersenneTwister(38)
        tot_amplo = rand(r4, 100:9999, 2000)
        p4 = rand(r4, Beta(8, 6), 2000)
        v_amplo = [rand(r4, Binomial(tot_amplo[i], p4[i])) for i in 1:2000]
        cal_amplo = calibration_check(v_amplo, tot_amplo;
                                      R = 120, tests = [:penultimate],
                                      rng = MersenneTwister(39))
        @test cal_amplo.verdict[1] !== :anticonservador
    end

    @testset "calibration_check — estrutura e Tables.jl" begin
        rng = seed("calibration_check — estrutura e Tables.jl")
        m = 400
        totals = rand(rng, 150:900, m)
        p = rand(rng, Beta(8, 6), m)
        votes = [rand(rng, Binomial(totals[i], p[i])) for i in 1:m]
        cal = calibration_check(votes, totals; R = 40, B = 99,
                                tests = [:last_digit, :rozenas], rng = rng)
        @test cal.tests == [:last_digit, :rozenas]
        @test cal.alpha == 0.05
        @test cal.m == m
        @test cal.R == 40
        @test all(0 .≤ cal.rejection_rate .≤ 1)
        @test all(cal.verdict .∈ Ref((:calibrado, :conservador, :anticonservador)))

        @test Tables.istable(typeof(cal))
        linhas = Tables.rows(cal)
        @test length(linhas) == 2
        @test keys(first(linhas)) == Tables.schema(cal).names
        @test Tables.columntable(cal).test == cal.tests
        @test Tables.columntable(cal).ratio ≈ cal.rejection_rate ./ cal.alpha

        # reprodutível com a mesma semente
        c1 = calibration_check(votes, totals; R = 40, B = 99,
                               tests = [:last_digit], rng = MersenneTwister(77))
        c2 = calibration_check(votes, totals; R = 40, B = 99,
                               tests = [:last_digit], rng = MersenneTwister(77))
        @test c1.rejection_rate == c2.rejection_rate

        io = IOBuffer(); show(io, MIME"text/plain"(), cal)
        @test occursin("Auto-diagnóstico", String(take!(io)))
    end

    @testset "_verdict — fronteiras" begin
        vd = ElectionForensics._verdict
        @test vd(0.05, 0.05, 0.01) === :calibrado
        @test vd(0.071, 0.05, 0.01) === :anticonservador   # > α + 2·se
        @test vd(0.069, 0.05, 0.01) === :calibrado
        @test vd(0.029, 0.05, 0.01) === :conservador       # < α − 2·se
        @test vd(0.031, 0.05, 0.01) === :calibrado
        @test vd(0.9, 0.05, 0.0) === :calibrado            # se = 0 ⇒ sem juízo
    end

    @testset "forensics_report — calibração integrada" begin
        rng = seed("forensics_report — calibração integrada")
        totals = rand(rng, 10:120, 1200)
        p = rand(rng, Beta(8, 6), 1200)
        votes = [rand(rng, Binomial(totals[i], p[i])) for i in 1:1200]
        buf = IOBuffer()
        forensics_report(votes, totals; B = 99, io = buf, rng = MersenneTwister(41))
        saida = String(take!(buf))
        @test occursin("Auto-diagnóstico", saida)
        @test occursin("anticonservador", saida)

        buf2 = IOBuffer()
        forensics_report(votes, totals; B = 99, calibrate = false, io = buf2,
                         rng = MersenneTwister(41))
        @test !occursin("Auto-diagnóstico", String(take!(buf2)))
    end

    @testset "_bb_logkernel — confere com a definição" begin
        # O núcleo omite log C(n,y); a diferença para o logpdf completo tem de
        # ser exatamente esse termo, para qualquer (α, β).
        for (y, n) in ((3, 10), (300, 500), (1, 2), (0, 7), (7, 7))
            for (a, b) in ((1.0, 1.0), (8.0, 6.0), (0.3, 0.4), (120.0, 35.0))
                k = ElectionForensics._bb_logkernel(y, n, a, b)
                completo = logpdf(BetaBinomial(n, a, b), y)
                @test completo - k ≈ first(logabsbinomial(n, y)) atol = 1e-9
            end
        end
        # α = β = 1 ⇒ a Beta-Binomial é uniforme em 0:n. Quem é constante é o
        # logpdf COMPLETO; o núcleo omite log C(n,y) e portanto não é.
        n = 12
        completos = [ElectionForensics._bb_logkernel(y, n, 1.0, 1.0) +
                     first(logabsbinomial(n, y)) for y in 0:n]
        @test all(≈(-log(n + 1)), completos)
        nucleos = [ElectionForensics._bb_logkernel(y, n, 1.0, 1.0) for y in 0:n]
        @test !all(≈(nucleos[1]), nucleos)
        @test nucleos ≈ reverse(nucleos)          # simétrico em α = β
    end

    @testset "_nelder_mead — valores analíticos" begin
        nm = ElectionForensics._nelder_mead
        # Rosenbrock: mínimo em (1,1)
        f(p) = (1 - p[1])^2 + 100(p[2] - p[1]^2)^2
        best, fv = nm(f, (-1.0, 1.0), -5.0, 5.0; maxiter = 3000, tol = 1e-14)
        @test best[1] ≈ 1.0 atol = 1e-2
        @test best[2] ≈ 1.0 atol = 2e-2
        @test fv < 1e-4
        # respeita a caixa
        g(p) = (p[1] - 10.0)^2 + (p[2] + 10.0)^2
        best2, _ = nm(g, (0.0, 0.0), -1.0, 1.0; maxiter = 500)
        @test -1.0 ≤ best2[1] ≤ 1.0
        @test -1.0 ≤ best2[2] ≤ 1.0
        @test best2[1] ≈ 1.0 atol = 1e-3     # empurrado contra a fronteira
    end

    @testset "fit_betabinomial_mixture — recupera mistura conhecida" begin
        r = MersenneTwister(101)
        m = 3000
        n = rand(r, 200:900, m)
        comp = rand(r, 1:2, m)
        # componente 1: média 0.400 (α=20,β=30); componente 2: média 0.750 (α=45,β=15)
        p = [comp[i] == 1 ? rand(r, Beta(20, 30)) : rand(r, Beta(45, 15)) for i in 1:m]
        y = [rand(r, Binomial(n[i], p[i])) for i in 1:m]

        mix = fit_betabinomial_mixture(y, n; rng = MersenneTwister(102))
        @test mix.components == 2
        medias = sort([mix.alpha[k] / (mix.alpha[k] + mix.beta[k]) for k in 1:2])
        @test medias[1] ≈ 0.400 atol = 0.02
        @test medias[2] ≈ 0.750 atol = 0.02
        @test sum(mix.weights) ≈ 1.0
        @test all(mix.weights .> 0)
        @test all(mix.alpha .> 0) && all(mix.beta .> 0)
        @test isfinite(mix.loglik) && isfinite(mix.bic)
        @test 1 ≤ mix.iterations ≤ 100
        @test mix.misfit ≥ 1.0                      # razão é sempre ≥ 1 por construção

        # BIC penaliza componentes a mais: com uma Beta única, escolhe L = 1
        r2 = MersenneTwister(103)
        n2 = rand(r2, 200:900, m)
        p2 = rand(r2, Beta(8, 6), m)
        y2 = [rand(r2, Binomial(n2[i], p2[i])) for i in 1:m]
        mix2 = fit_betabinomial_mixture(y2, n2; rng = MersenneTwister(104))
        @test mix2.components == 1
        @test mix2.alpha[1] / (mix2.alpha[1] + mix2.beta[1]) ≈ 8/14 atol = 0.02

        @test_throws ArgumentError fit_betabinomial_mixture(y, n; max_components = 0)
        @test_throws ArgumentError fit_betabinomial_mixture(y, n; log_bounds = (5.0, 1.0))
    end

    @testset "rozenas_test — nula Beta-Binomial" begin
        rng = seed("rozenas_test — nula Beta-Binomial")
        m = 1200
        totals = rand(rng, 150:900, m)
        p = rand(rng, Beta(8, 6), m)
        votes = [rand(rng, Binomial(totals[i], p[i])) for i in 1:m]

        rk = rozenas_test(votes, totals; B = 199, null = :kernel,
                          rng = MersenneTwister(51))
        rb = rozenas_test(votes, totals; B = 199, null = :betabinomial,
                          rng = MersenneTwister(51))
        @test rk.null === :kernel
        @test rb.null === :betabinomial
        @test rk.mixture === nothing
        @test rb.mixture isa BetaBinomialMixture
        @test rb.total_observed == rk.total_observed     # a estatística não muda
        @test 0 < rb.pvalue ≤ 1
        @test length(rb.qvalues) == length(rb.fractions)

        # reprodutível
        rb2 = rozenas_test(votes, totals; B = 199, null = :betabinomial,
                           rng = MersenneTwister(51))
        @test rb.pvalue == rb2.pvalue

        @test_throws ArgumentError rozenas_test(votes, totals; null = :bogus)

        io = IOBuffer(); show(io, MIME"text/plain"(), rb)
        @test occursin("Beta-Binomial", String(take!(io)))
    end

    @testset "Aqua" begin
        Aqua.test_all(ElectionForensics)
    end
end
