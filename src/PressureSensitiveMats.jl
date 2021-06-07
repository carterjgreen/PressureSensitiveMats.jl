module PressureSensitiveMats

using Distances, DSP, FFTW, InvertedIndices, Peaks, StatsBase, Statistics

export est_br, est_br_fft, est_br_fft2
export PCC, PCC2, SNR_MAX, EGC, MRC_PSD, get_weights, combiner
export Holtz, Solei, move_detect
export breath_availability, occupancy_detection
export active_sensors,
       apply2seg,
       choose_ref,
       estimate_snr,
       extract_ref,
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