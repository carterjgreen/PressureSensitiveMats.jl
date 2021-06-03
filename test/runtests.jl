using Base: Number
using PressureSensitiveMats
using Test

@testset "PressureSensitiveMats.jl" begin
    @test sfm(zeros(512)) isa Number

end

@testset "Combiners" begin
    multi = fill(313.0, (35000, 72)) .+ randn((35000, 72))
    @test combiner(multi) isa AbstractVector
end

@testset "Occupancy Detection" begin
    @test occupancy_detection(fill(313.0, (35000, 72)) .+ randn((35000, 72))) == zeros(35000)
end

@testset "Motion Detection" begin
    single = fill(313.0, 35000) .+ randn(35000)
    multi = fill(313.0, (35000, 72)) .+ randn((35000, 72))
    @test move_detect(single) == zeros(35000)
    @test move_detect(multi) == zeros(35000)
    @test move_detect(single, Solei(), height=50, weight=50) == zeros(35000)
end