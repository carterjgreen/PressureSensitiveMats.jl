struct Holtz end
struct Solei end

"""
    interval_merging(onset::BitVector, offset::BitVector)

Given two BitVectors signalling offsets and onsets, this function merges the two 
to create one indicator signal.
"""
function interval_merging(onset::BitVector, offset::BitVector)
    out = falses(length(onset))
    if any(onset)
        for i in findall(onset)
            ind = findnext(offset, i)
            if !isnothing(ind)
                out[i:ind] .= true
            end
        end
    end
    return out
end

"""
    window_sizes(x::AbstractVector, w_max=300, thresh=100)

Function to calculate window sizes based on location of peaks in a combined PSM signal.

# Arguments:
 - `w_max=300` : Maximum window size.
 - `thresh=100` : Minimum peak magnitude threshold.
"""
function window_sizes(x::AbstractVector{<:Number}, w_max::Integer = 300,
                      thresh::Number = 100)
    windows = fill(w_max, length(x))
    inds = argmaxima(x, 20)
    filter!(i -> x[i] > thresh, inds)
    diffs = diff(inds)

    for i in 2:(length(inds) - 1)
        val = min(diffs[i - 1], w_max)
        new_w = @view windows[inds[i]:inds[i + 1]]
        fill!(new_w, val)
    end
    return windows
end

"""
    move_detect(::Holtz, x::AbstractVector; kwargs...)

Implements movement detection from Identifying Movement Onset Times for a Bed-Based 
Pressure Sensor Array - 2006, Holtzman.

# Arguments:
 - `L` : Length of window for moving moving_stats.
 - `κ` : Magnitude of the control signal.
 - `min_samples` : The minimum number of samples that a movement is allowed to be.

"""
function move_detect(::Holtz, x::AbstractVector{<:Number}; L::Integer = 300, κ::Real = 3,
                     min_samples::Integer = 3)
    # Movement detection from Holtzman, 2006
    onset = falses(length(x))
    offset = falses(length(x))
    # Handle empty beds

    moving_avg, moving_var = moving_stats(x, L)
    mvar_diff = diff(moving_var) # diff([0; mvar]) optionally
    control = κ .* sqrt.(moving_var)

    ucl = moving_avg .+ control
    lcl = moving_avg .- control

    @views for i in 300:(length(x) - (min_samples + 1))
        val = x[i:(i + min_samples - 1)]
        if all((val .< lcl[i:(i + min_samples - 1)]) .|
               (val .> ucl[i:(i + min_samples - 1)]))
            onset[i] = true
        end
    end

    # Flag offsets
    for i in findall(onset) # location of trues
        ind = findfirst(x -> x < 0.1, @view mvar_diff[i:end])
        !isnothing(ind) ? offset[i + ind] = true : nothing
    end

    return interval_merging(onset, offset)
end

function move_detect(x::AbstractMatrix{<:Number}; L::Integer = 300, κ::Real = 3,
                     min_samples::Integer = 3)
    out = falses(size(x))
    if mean(occupancy_detection(x)) > 0.75
        for (i, col) in enumerate(eachcol(x))
            # Apply detection on each column/sensor
            out[:, i] = move_detect(Holtz(), col, L = L, κ = κ, min_samples = min_samples)
        end
    end
    return vec(sum(out, dims = 2)) .>= 2
end

"""
    move_detect(::Solei, x::AbstractVector; kwargs...)

Implements movement detection from Movement Detection with Adaptive Window Length for 
Unobtrusive Bed-based Pressure-Sensor Array - 2017, Soleimani

# Arguments:
 - `L` : Length of window for moving moving_stats.
 - `α=-0.029` : Threshold of first differences for the offset
 - `κ=3` : Magnitude of the control signal.
 - `min_samples=2` : The minimum number of samples that a movement is allowed to be.
 - `height=50` : Subject height in cm. Default to cancel out weight in ρ calculation.
 - `weight=50` : Subject weight in kg. Default to cancel out height in ρ calculation.

"""
function move_detect(::Solei,
                     x::AbstractVector{<:Number};
                     α::Real = -0.029, κ::Integer = 3, min_samples::Integer = 2,
                     height::Real = 50, weight::Real = 50)
    # Movement detection from Soleimani, 2017
    # Expects reference sensor that is band-pass filtered
    # Set height and weight to 50 to ignore the last threshold
    onset = falses(length(x))
    offset = falses(length(x))

    windows = window_sizes(x)
    ρ = 3 * (weight + height - 100)
    moving_avg, moving_var = moving_stats(x, windows)
    control = κ .* sqrt.(moving_var)

    ucl = moving_avg .+ control
    ucl_diff = diff([0; ucl])
    lcl = moving_avg .- control

    @views for i in 300:(length(x) - (min_samples + 1))
        val = x[i:(i + min_samples - 1)]
        if all((val .< lcl[i:(i + min_samples - 1)]) .|
               (val .> ucl[i:(i + min_samples - 1)]))
            onset[i] = true
        end
    end

    # Flag offsets
    @views for i in findall(onset) # location of trues
        for j in i:(length(onset) - 2)
            if all(ucl_diff[j:(j + 2)] .< α)
                offset[j] = true
                break
            end
        end
    end

    indicator = interval_merging(onset, offset)
    return indicator .& (moving_var .> ρ)
end

move_detect(x; kwargs...) = move_detect(Holtz(), x; kwargs...)
