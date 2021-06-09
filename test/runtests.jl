using Base: Number
using PressureSensitiveMats
using Test

@testset "Utils" begin
    ma1, mv1 = moving_stats(ones(5000), 3)
    ma2, mv2 = moving_stats(ones(5000), fill(3, 5000))
    truth = [
        21  10  24   5  19   8  14   3
        9  23  12  18   7  13   2  16
       22  11  17   6  20   1  15   4
       45  34  48  29  43  32  38  27
       33  47  36  42  31  37  26  40
       46  35  41  30  44  25  39  28
       69  58  72  53  67  56  62  51
       57  71  60  66  55  61  50  64
       70  59  65  54  68  49  63  52]

    @test ma1[3:end] ≈ ones(4998)
    @test mv1[5:end] ≈ zeros(4996)

    @test ma2 ≈ ones(5000)
    @test mv2 ≈ zeros(5000)

    @test mat_shape(ones(72)) isa AbstractVector
    @test mat_shape(ones(5000, 72)) isa AbstractMatrix

    @test sfm(zeros(512)) isa Number
    @test active_sfm(randn(35000, 72), 300) isa AbstractVector

    @test reshape_psm(ones(Int, 400, 72) .* (1:72)')[:, :, rand(1:72)] == truth
end

@testset "Combiners" begin
    multi = fill(313.0, (35000, 72)) .+ randn((35000, 72))

    @test combiner(SNR_MAX(), multi) isa AbstractVector
    @test combiner(PCC(), multi) isa AbstractVector
    @test combiner(PCC2(), multi) isa AbstractVector
    @test combiner(EGC(), multi) isa AbstractVector
    @test apply2seg(s -> combiner(MRC_PSD(10), s), multi, 300) isa AbstractVector
end

@testset "Occupancy Detection" begin
    multi = fill(313.0, (35000, 72)) .+ randn((35000, 72))

    multi_occ = occupancy_detection(multi)
    multi_occ_dist = occupancy_detection(multi, max_dist=true)

    @test multi_occ isa BitVector
    @test multi_occ == zeros(35000)

    @test multi_occ_dist isa BitVector
    @test multi_occ_dist == zeros(35000)
end

@testset "Motion Detection" begin
    single = fill(313.0, 35000) .+ randn(35000)
    multi = fill(313.0, (35000, 72)) .+ randn((35000, 72))

    md_s = move_detect(single, min_samples=10)
    md_m = move_detect(multi, min_samples=10)
    md_solei = move_detect(Solei(), single, min_samples=10)

    @test md_s isa BitVector
    @test md_s == zeros(35000)

    @test  md_m isa BitVector
    @test md_m == zeros(35000)

    @test md_solei isa BitVector
    @test md_solei == zeros(35000)
end