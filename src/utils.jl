# Utility functions

ord = [
	21, 17, 13, 
	9, 5, 1, 
	22, 18, 14,
	10, 6, 2,
	23, 19, 15,
	11, 7, 3, 
	24, 20, 16, 
	12, 8, 4]

function estimate_snr(x::AbstractVector; fs=10)
    # Estimate SNR for a 30-102.4s segment
    pow = power(periodogram(x, nfft=1024, window=hanning))
    ps = argmax(pow)
    n = mean(pow[Not(ps, 2ps)]) * 0.73
    w = (pow[ps] * fs / 1024 - n) / n
    return w
end

"""
    moving_stats(x::AbstractVector, L::Int)

Return moving average and variance of x using a window length L
"""
function moving_stats(x::AbstractVector{T}, L::Int) where T
    # Use filters to calculate running mean and var
    avg_filt = ones(T, L) ./ L
    var_filt = ones(T, L) ./ (L - 1) # unbiased
    
    moving_avg = filt(avg_filt, 1, x)
    moving_var = filt(var_filt, 1, (x .- moving_avg) .^ 2)
    
    return moving_avg, moving_var
end

"""
    moving_stats(x::AbstractVector, w::AbstractVector{Int})

Return moving average and variance of x using a window lengths w
"""
function moving_stats(x::AbstractVector{T}, w::AbstractVector{Int}) where T
    moving_avg = zeros(T, length(x))
    moving_var = zeros(T, length(x))
    @views for i in eachindex(x)
        win = i < w[i] ? i : w[i]
        moving_avg[i] = mean(x[i-win+1:i])
        moving_var[i] = sum((x[i-win+1:i] .- moving_avg[i-win+1:i]) .^ 2) / (w[i]-1)
    end
    return moving_avg, moving_var
end

"""
    apply2seg(f::Function, x::AbstractMatrix{T} , n::Integer)

Convenience function to apply a function f to segments of length n to the matrix x.
"""
function apply2seg(f::Function, x::AbstractMatrix{T} , n::Integer) where T
    # Assumes a mat comes in and a vec goes out
    # Not much better than a mapreduce but it includes last segment
    ra = 1:n:size(x, 1)-n
    out = Vector{T}(undef, size(x, 1))
    for i in ra
        out[i:i+n-1] = f(@view x[i:i+n-1, :])
    end
    out[ra.stop+n:end] = f(@view x[ra.stop+n:end, :])
    return out
end

"""
    mat_shape(x, n)

Reorders the columns of a matrix to be in the proper sensor order.
"""
function mat_shape(x::AbstractVector, n=3)
    return view(x, mapreduce(i -> ord .+ (24 * i), vcat, 0:n-1))
end

function mat_shape(x::AbstractMatrix, n=3)
    return view(x, :, mapreduce(i -> ord .+ (24 * i), vcat, 0:n-1))
end

"""
	reshape_psm(psm, n_mats)
Takes in a PSM and reshapes it to be in the correct order.

ex. Nx72 -> 9x8xN 

ex. Nx24 -> 3x8xN
"""
function reshape_psm(x::AbstractMatrix, n=div(size(x, 2), 24))
	new_shape = reshape(mat_shape(x, n)', 3, 8, n, :) 
	return reduce(vcat, eachslice(new_shape, dims=3))
end

"""
    active_sensors(x, thresh)

    Return the indices of sensors whose mean value is greater than thresh.
"""
active_sensors(x::AbstractMatrix, thresh=0.4*2046) = vec(mean(x, dims=1) .> thresh)

"""
    choose_ref(x)

Return the index of the sensor with the greatest power.
"""
function choose_ref(x::AbstractMatrix)
    # Choose 
    out = sum(abs2, x, dims=1) |> argmax
    return out[2]
end

extract_ref(x::AbstractMatrix) = x[:, choose_ref(x)]

"""
    polarity_flip(x)

Flips the polarity of sensors based on their PCC with the reference sensor.
"""
polarity_flip(x) = sign.(cor(x, @view x[:, choose_ref(x)]))' .* x

"""
    sfm(x)
Returns the spectral flatness measure of a signal.
"""
function sfm(x)
    s = abs2.(fft(x))
    return geomean(s) / mean(s)
end

"""
    active_sfm(x, n, thresh)

Return the indices of sensors with a mean spectral flatness measure, taken in segments of
length n, above the threshold in dB.
"""
function active_sfm(x::AbstractMatrix, n, thresh=-50)
    fo(sig) = mapreduce(i -> sfm(sig[i:i+n-1]), vcat, 1:n:size(x, 1)-n)
    s = mapreduce(fo, hcat, eachcol(x))
    actives = pow2db.(s) .> thresh
    return vec(mean(actives, dims=1) .> 0.4)
end
