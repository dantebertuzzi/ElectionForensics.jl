using ElectionForensics
using Test
using Random
using Distributions

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
        @test penultimate_digit(1234) == 3
        @test penultimate_digit(10) == 1
        @test_throws ArgumentError penultimate_digit(5)
    end

    @testset "coarse_fractions" begin
        fr = coarse_fractions(4)
        @test fr == sort([1//2, 1//3, 2//3, 1//4, 3//4])
        @test all(f -> 0 < f < 1, coarse_fractions(10))
        @test allunique(coarse_fractions(10))
        @test_throws ArgumentError coarse_fractions(1)
    end

    rng = MersenneTwister(2026)

    @testset "Benford — dados log-uniformes conformes" begin
        # log-uniforme em [1, 10^5) segue Benford exatamente
        x = [round(Int, 10.0^(5rand(rng))) for _ in 1:5000]
        x = [xi for xi in x if xi ≥ 10]
        r1 = benford_test(x; digit = 1)
        r2 = benford_test(x; digit = 2)
        @test r1.pvalue > 0.001
        @test r2.pvalue > 0.001
        @test r1.conformity in (:conforme, :aceitavel)
        @test sum(r1.counts) == r1.n
        @test isapprox(sum(r1.expected), 1.0; atol = 1e-10)
        @test isapprox(sum(r2.expected), 1.0; atol = 1e-10)
    end

    @testset "Benford — dados fabricados não conformes" begin
        # primeiro dígito uniforme em 1:9 viola fortemente a 1BL
        x = [d * 10^rand(rng, 1:4) + rand(rng, 0:9) for d in rand(rng, 1:9, 3000)]
        r = benford_test(x; digit = 1)
        @test r.pvalue < 1e-6
        @test r.conformity == :nao_conforme
    end

    @testset "último dígito — limpo vs. fabricado" begin
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

    @testset "Rozenas — validação de entradas" begin
        @test_throws DimensionMismatch rozenas_test([1, 2], [10])
        @test_throws ArgumentError rozenas_test([5], [0])
        @test_throws ArgumentError rozenas_test([11], [10])
        @test_throws ArgumentError rozenas_test(Int[], Int[])
    end

    @testset "forensics_report" begin
        m = 800
        totals = rand(rng, 100:600, m)
        votes = [rand(rng, Binomial(totals[i], 0.55)) for i in 1:m]
        res = forensics_report(votes, totals; B = 99, rng = rng,
                               io = IOBuffer())
        @test res.benford isa BenfordResult
        @test res.last_digit isa LastDigitResult
        @test res.rozenas isa RozenasResult
    end
end
