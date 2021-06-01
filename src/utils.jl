# Utility functions

ord = (
    21, 9, 22, 10, 23, 11, 24, 12, 
    17, 5, 18, 6, 19, 7, 20, 8, 
    13, 1, 14, 2, 15, 3, 16, 4)

function moving_stats(x::AbstractVector{T}, L::Int) where T
    # Use filters to calculate running mean and var
    avg_filt = ones(T, L) ./ L
    var_filt = ones(T, L) ./ (L - 1) # unbiased
    
    moving_avg = filt(avg_filt, 1, x)
    moving_var = filt(var_filt, 1, (x .- moving_avg) .^ 2)
    
    return moving_avg, moving_var
end

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

function apply2seg(x::AbstractMatrix{T}, f::Function, n::Integer) where T
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

function mat_shape(x::AbstractVector, n=3)
    return view(x, mapreduce(i -> ord .* i, vcat, 1:n))
end

function mat_shape(x::AbstractMatrix, n=3)
    return view(x, :, mapreduce(i -> ord .* i, vcat, 1:n))
end

active_sensors(x::AbstractMatrix, thresh=0.4*2046) = vec(mean(x, dims=1) .> thresh)

function choose_ref(x::AbstractMatrix)
    # Choose 
    out = sum(abs2, x, dims=1) |> argmax
    return out[2]
end

extract_ref(x::AbstractMatrix) = x[:, choose_ref(x)]

polarity_flip(x) = sign.(cor(x, @view x[:, choose_ref(x)]))' .* x

function sfm(x)
    s = abs2.(fft(x))
    return geomean(s) / mean(s)
end

"asfsa"