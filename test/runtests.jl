using Base: Number
using PressureSensitiveMats
using Test

@testset "PressureSensitiveMats.jl" begin
    @test sfm(zeros(512)) isa Number
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