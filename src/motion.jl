struct Holtz end
struct Solei end

function interval_merging(onset::BitVector, offset::BitVector)
    out = falses(length(onset))
    for i in findall(onset)
        ind = findnext(offset, i)
        out[i:ind] .= true
    end
    return out
end

function window_sizes(x::AbstractVector, w_max=300, thresh=100)
    windows = fill!(Vector{Int}(undef, length(x)), w_max)
    inds = argmaxima(x, 20)
    filter!(i -> x[i] > thresh, inds)
    diffs = diff(inds)
    
    for i in 2:length(inds)-1
        val = min(diffs[i-1], w_max)
        new_w = @view windows[inds[i]:inds[i+1]]
        fill!(new_w, val)
    end
    return windows
end

function move_detect(x::AbstractVector, ::Holtz; L=300, κ=3, min_samples=3)
    # Movement detection from Holtzman, 2006
    moving_avg, moving_var = moving_stats(x, L)
    mvar_diff = diff(moving_var) # diff([0; mvar]) optionally
    control = κ .* sqrt.(moving_var)
    
    ucl = moving_avg .+ control
    lcl = moving_avg .- control
    
    onset = falses(length(x))
    offset = falses(length(x))

    @views for i in 300:length(x)-(min_samples+1)
        val = x[i:i+min_samples-1]
        if all((val .< lcl[i:i+min_samples-1]) .| (val .> ucl[i:i+min_samples-1]))
            onset[i] = true
        end
    end
    
    # Flag offsets
    for i in findall(onset) # location of trues
        ind = findfirst(x -> x < 0.1, @view mvar_diff[i:end])
        offset[i + ind] = true
    end
    
    return interval_merging(onset, offset)
end

function move_detect(x::AbstractMatrix; L=300, κ=3, min_samples=3)
    out = falses(size(x))
    for (i, col) in enumerate(eachcol(x))
    # Apply detection on each column/sensor
        out[:, i] = move_detect(col, Holtz(), L=L, κ=κ, min_samples=min_samples)
    end
    
    return vec(sum(out, dims=2)) .>= 2
end

function move_detect(x::AbstractVector, ::Solei; α=-0.029, κ=3, min_samples=2, height=183, weight=93)
    # Movement detection from Soleimani, 2017
    # Expects reference sensor that is band-pass filtered
    windows = window_sizes(x)
    ρ = 3 * (weight + height - 100)
    moving_avg, moving_var = moving_stats(x, windows)
    control = κ .* sqrt.(moving_var)
    
    ucl = moving_avg .+ control
    ucl_diff = diff([0; ucl])
    lcl = moving_avg .- control
    
    onset = falses(length(x))
    offset = falses(length(x))

    @views for i in 300:length(x)-(min_samples+1)
        val = x[i:i+min_samples-1]
        if all((val .< lcl[i:i+min_samples-1]) .| (val .> ucl[i:i+min_samples-1]))
            onset[i] = true
        end
    end
    
#     Flag offsets
    @views for i in findall(onset) # location of trues
        for j in i:length(onset)-2
            if all(ucl_diff[j:j+2] .< α)
                offset[j] = true
                break
            end
        end
    end
    
    indicator = interval_merging(onset, offset)
    return indicator .& (moving_var .> ρ)
end

move_detect(x, kwargs...) = move_detect(x, Holtz(), kwargs...)
