using BenchmarkTools, PressureSensitiveMats, Random

SUITE = BenchmarkGroup()

SUITE["motion"] = BenchmarkGroup()
SUITE["combiners"] = BenchmarkGroup()
SUITE["occupancy"] = BenchmarkGroup()

# Inits
single = fill(313.0, 35000) .+ randn(35000)
multi = fill(313.0, (35000, 72)) .+ randn((35000, 72))
n = 10 * 60 * 60 * 48 # More than 48 hours crashes BenchmarkTools.jl
two_days = fill(312.0, (n, 72)) .+ randn((n, 72))

# Motion
SUITE["motion"]["Single Holtzmann"] = @benchmarkable move_detect($single, min_samples=10)
SUITE["motion"]["Multi Holtzmann"] = @benchmarkable move_detect($multi, min_samples=10)
SUITE["motion"]["Soleimani"] = @benchmarkable move_detect(Solei(), $single, min_samples=10)

# Occupancy
SUITE["occupancy"]["regular"] = @benchmarkable occupancy_detection($multi)
SUITE["occupancy"]["dist method"] = @benchmarkable occupancy_detection($multi, max_dist=true)

SUITE["occupancy"]["two days regular"] = @benchmarkable occupancy_detection($two_days)
SUITE["occupancy"]["two days dist method"] = @benchmarkable occupancy_detection($two_days, max_dist=true)
# Combiners
SUITE["combiners"]["SNR-MAX"] = @benchmarkable combiner(SNR_MAX(), $multi)
SUITE["combiners"]["Pearson Correlation Coefficient"] = @benchmarkable combiner(PCC(), $multi)
SUITE["combiners"]["Pearson Correlation Coefficient v2"] = @benchmarkable combiner(PCC2(), $multi)
SUITE["combiners"]["Equal Gain Combining"] = @benchmarkable combiner(EGC(), $multi)