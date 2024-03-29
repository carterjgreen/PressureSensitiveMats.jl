using Base: Number
using PressureSensitiveMats
using Aqua, Test

# TODO Add type tests for float32
@testset "Utils" begin
    ma1, mv1 = moving_stats(ones(5000), 3)
    ma2, mv2 = moving_stats(ones(5000), fill(3, 5000))
    truth = [21 9 22 10 23 11 24 12
        17 5 18 6 19 7 20 8
        13 1 14 2 15 3 16 4
        45 33 46 34 47 35 48 36
        41 29 42 30 43 31 44 32
        37 25 38 26 39 27 40 28
        69 57 70 58 71 59 72 60
        65 53 66 54 67 55 68 56
        61 49 62 50 63 51 64 52]

    active = hcat(fill(1000, 40), zeros(40))

    @test ma1[3:end] ≈ ones(4998)
    @test mv1[5:end] ≈ zeros(4996)

    @test ma2 ≈ ones(5000)
    @test mv2 ≈ zeros(5000)

    @test mat_shape(ones(72)) isa AbstractVector{<:Number}
    @test mat_shape(ones(5000, 72)) isa AbstractMatrix{<:Number}

    @test sfm(zeros(512)) isa Number
    @test active_sfm(randn(35000, 72), 300) isa AbstractVector
    @test active_sensors(active) == Bool[1, 0]
    @test extract_ref(active) == active[:, 1]

    @test reshape_psm(ones(Int, 400, 72) .* (1:72)')[:, :, rand(1:72)] == truth
end

@testset "Breathing Rate" begin
    breath = cos.(2π * 12 / 60 / 10 .* collect(1000:4001))

    @test 11.5 < est_br(breath) < 12.5
    @test 11.8 < est_br_fft(breath) < 12.2
    @test 11.8 < est_br_fft2(breath) < 12.2
end

@testset "Combiners" begin
    multi = fill(313.0, (35000, 72)) .+ randn((35000, 72))

    @test combiner(SNR_MAX(), multi) isa AbstractVector{<:Number}
    @test snrmax(multi) == combiner(SNR_MAX(), multi)
    @test combiner(multi) == combiner(SNR_MAX(), multi)

    @test combiner(PCC(), multi) isa AbstractVector{<:Number}
    @test combiner(PCC2(), multi) isa AbstractVector{<:Number}
    @test pcc(multi) == combiner(PCC2(), multi)

    @test combiner(EGC(), multi) isa AbstractVector{<:Number}
    @test egc(multi) == combiner(EGC(), multi)

    @test selection(multi) isa AbstractVector{<:Number}

    @test combiner(SNR_MAX(), multi, 5) isa AbstractVector{<:Number}
    @test apply2seg(s -> combiner(MRC_PSD(), s), multi, 300) isa AbstractVector{<:Number}
    @test apply2seg(s -> combiner(MRC_PSD(10), s), multi, 300) isa AbstractVector{<:Number}
    @test apply2seg(s -> mrc(s, 10), multi, 300) isa AbstractVector{<:Number}
end

@testset "Occupancy Detection" begin
    multi = fill(313.0, (35000, 72)) .+ randn((35000, 72))
    multi2 = copy(multi)
    multi2[1000:4000, :] .= 2000.0

    multi_occ = occupancy_detection(multi)
    multi_occ2 = occupancy_detection(multi2)
    multi_occ_dist = occupancy_detection(multi; max_dist = true)
    multi_occ_dist2 = occupancy_detection(multi2; max_dist = true)

    @test multi_occ isa BitVector
    @test multi_occ == zeros(35000)
    @test sum(multi_occ2) == 3001

    @test multi_occ_dist isa BitVector
    @test multi_occ_dist == zeros(35000)
    @test multi_occ_dist2 != multi_occ2 # can tell that it's not real
end

@testset "Motion Detection" begin
    single = fill(313.0, 35000) .+ randn(35000)
    multi = fill(313.0, (35000, 72)) .+ randn((35000, 72))

    single2 = copy(single)
    single2[1000:1200] .+= 1000

    md_s = move_detect(single; min_samples = 10)
    md_s2 = move_detect(single2; min_samples = 10)
    md_m = move_detect(multi; min_samples = 10)
    md_solei = move_detect(Solei(), single; min_samples = 10)

    @test md_s isa BitVector
    @test md_s == zeros(35000)
    @test sum(md_s2) == 301

    @test md_m isa BitVector
    @test md_m == zeros(35000)

    @test md_solei isa BitVector
    @test md_solei == zeros(35000)
end

@testset "Code Quality" begin
    Aqua.test_ambiguities(PressureSensitiveMats; recursive = false)
    # Aqua.test_all includes Base and Core in ambiguity testing
    Aqua.test_all(PressureSensitiveMats; ambiguities = false)
end
