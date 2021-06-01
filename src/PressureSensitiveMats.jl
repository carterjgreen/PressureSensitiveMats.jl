module PressureSensitiveMats

using Distances, DSP, FFTW, Peaks, StatsBase, Statistics

include("breathing.jl")
include("combiners.jl")
include("motion.jl")
include("occupancy.jl")
include("utils.jl")

end
