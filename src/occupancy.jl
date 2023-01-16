"""
    min_occupancy(x, m=300)

Given an indicator signal this function uses run-length encoding to enforce minimum
lengths. If an indicator segment is less than the minimum it is added to the length
of the previous segment. 
"""
function min_occupancy(x::AbstractVector{<:Number}, m::Integer=300)
    # Use run-length encoding to find lengths and enforce minimum
    a, b = rle(x)
    a_new, b_new = Bool[], Int[]
    for i in eachindex(b)
        if b[i] > m || isempty(b_new)
            push!(a_new, a[i])
            push!(b_new, b[i])
        else
            b_new[end] += b[i]
        end
    end
    return BitVector(inverse_rle(a_new, b_new))
end

"""
    breath_availability(x)

Implements the breathing availbaility indicator from Long-Term Sleep Assessment by 
Unobtrusive Pressure Sensor Arrays - 2018, Soleimani. 

Expects an array with dimensions n_mats x 8 x samples. Returns the maximum distance 
between mattress sensors and the mode of those distances.
"""
function breath_availability(x::AbstractArray{T, 3}) where T<:Number
    # Needs unfiltered data in a 9x8xN format in the proper order
    out = Vector{T}(undef, size(x, 3))
    for i in eachindex(out)
        active = findall(view(x, :, :, i) .> (2046 * 0.4))
        if !isempty(active)
            out[i] = pairwise(Euclidean(), Tuple.(active)) |> maximum
        else
            out[i] = Inf
        end
    end
    return mode(out), out
end

"""
    occupancy_detection(x::AbstractVector; kwargs)

Implements signle sensor simple occupancy detection from Long-Term Sleep Assessment by 
Unobtrusive Pressure Sensor Arrays - 2018, Soleimani.

# Arguments
 - `β=350` : Average value of a sensor from an empty bed.
 - `m=300` : Minimum number of samples to be declared an occupancy.
 - `n=4` : Sensitivity parameter
"""
function occupancy_detection(x::AbstractVector{<:Number}; β::Real=350, m::Integer=300, n::Real=4)
    τ = cityblock(extrema(x)...) / n # Sareh uses n=4 in her code
    occupancy = (x .> (β + τ))
    return any(occupancy) ? min_occupancy(occupancy, m) : occupancy
end

"""
    occupancy_detection(x::AbstractMatrix; kwargs)

Implements simple occupancy detection from Long-Term Sleep Assessment by Unobtrusive
Pressure Sensor Arrays - 2018, Soleimani.

# Arguments
 - `max_dist=false` : Use breathing availability signal from the maximum distance method.
 - `β=25000` : Average value of the sum of all sensors from an empty bed.
 - `m=300` : Minimum number of samples to be declared an occupancy.
 - `n=4` : Sensitivity parameter
"""
function occupancy_detection(x::AbstractMatrix{<:Number}; max_dist::Bool=false, β::real=25000, m::Integer=300, n::Real=4)
    Z = sum(x, dims=2) |> vec
    τ = cityblock(extrema(Z)...) / n # Sareh uses n=4 in her code
    
    if max_dist # Sareh's post-processing
        k, breath = breath_availability(reshape_psm(x))
        occupancy = (Z .> (β + τ)) .& (0.85*k .< breath .< 1.15*k)
    else
        occupancy = (Z .> (β + τ))
    end
    return any(occupancy) ? min_occupancy(occupancy, m) : occupancy
end
