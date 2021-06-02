module PressureSensitiveMats

using Distances, DSP, FFTW, InvertedIndices, Peaks, StatsBase, Statistics

export est_br, est_br_fft, est_br_fft2
export PCC, PCC2, SNR_MAX, EGC, MRC_PSD, get_weights, combiner
export move_detect
export active_sensors,
       apply2seg,
       choose_ref,
       estimate_snr,
       mat_shape, 
       moving_stats,
       polarity_flip,
       sfm
       
include("breathing.jl")
include("combiners.jl")
include("motion.jl")
include("occupancy.jl")
include("utils.jl")

end
