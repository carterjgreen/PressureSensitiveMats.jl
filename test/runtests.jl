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

    @test occupancy_detection(multi) == zeros(35000)
end

@testset "Motion Detection" begin
    single = fill(313.0, 35000) .+ randn(35000)
    multi = fill(313.0, (35000, 72)) .+ randn((35000, 72))

    @test move_detect(single, min_samples=10) == zeros(35000)
    @test move_detect(multi, min_samples=10) == zeros(35000)
    @test move_detect(single, Solei(), min_samples=10) == zeros(35000)
end