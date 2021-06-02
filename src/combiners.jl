abstract type Combiner end
abstract type CorPol <: Combiner end

struct PCC <: Combiner end
struct PCC2 <: Combiner end
struct SNR_MAX <: Combiner end
struct EGC <: CorPol end

struct MRC_PSD <: CorPol 
    fs:Int
end

MRC_PSD(fs=10) = MRC_PSD(fs)

"""
    get_weights(x, ref, comb)

Returns the weights for a section of the PSM by the method comb
"""
function get_weights(x::AbstractMatrix, ref::Int, comb::PCC)
    return cor(x, view(x, :, ref)) |> vec
end

function get_weights(x::AbstractMatrix, ref::Int, comb::PCC2)
    pows = sum(abs2, x, dims=1) .|> sqrt
    cors = @views sum(x .* x[:, ref], dims=1)
    w = cors ./ (pows .* pows[ref]) |> vec
    return w
end

function get_weights(x::AbstractMatrix, ref::Int, comb::SNR_MAX)
    cors = @views sum(x .* x[:, ref], dims=1)
    w = cors ./ cors[ref] |> vec
    return w
end

get_weights(x::AbstractMatrix, comb::Combiner) = get_weights(x, choose_ref(x), comb)

function get_weights(x::AbstractMatrix{T}, comb::MRC_PSD) where T
    w = Vector{T}(undef, size(x, 2))
    for (i, s) in enumerate(eachcol(x))
        w[i] = estimate_snr(s, fs=comb.fs)
    end
    return w
end

function get_weights(x::AbstractMatrix{T}, comb::EGC) where T
    return ones(T, size(x, 2))
end

"""
    combiner(x, comb)

Gets weights from comb and applies them to x
"""
function combiner(x::AbstractMatrix, comb::Combiner)
    w = get_weights(x, comb)
    return x * w
end

function combiner(x::AbstractMatrix, comb::CorPol)
    w = get_weights(x, comb)
    return polarity_flip(x) * w
end

combiner(x::AbstractMatrix) = combiner(x, SNR_MAX())