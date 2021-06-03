using Base: Number
using PressureSensitiveMats
using Test

@testset "PressureSensitiveMats.jl" begin
    @test sfm(zeros(512)) isa Number
end
