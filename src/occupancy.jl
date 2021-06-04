function min_occupancy(x, m=300)
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

function breath_availability(x::AbstractArray{T, 3}) where T
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

function occupancy_detection(x::AbstractVector; β=312, m=300, n=4)
    τ = cityblock(extrema(x)...) / n # Sareh uses n=4 in her code
    occupancy = (x .> (β + τ))
    return any(occupancy) ? min_occupancy(occupancy, m) : occupancy
end

function occupancy_detection(x::AbstractMatrix; max_dist=false, β=22555, m=300, n=4)
    Z = sum(x, dims=2) |> vec
    τ = cityblock(extrema(Z)...) / n # Sareh uses n=4 in her code
    
    if max_dist # Sareh's post-processing
        k, breath = breath_availability(reshape(mat_shape(x)', 9, 8, :))
        occupancy = (Z .> (β + τ)) .& (0.85*k .< breath .< 1.15*k)
    else
        occupancy = (Z .> (β + τ))
    end
    return any(occupancy) ? min_occupancy(occupancy, m) : occupancy
end
