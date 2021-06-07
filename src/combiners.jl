abstract type Combiner end
abstract type CorPol <: Combiner end

struct PCC <: Combiner end
struct PCC2 <: Combiner end
struct SNR_MAX <: Combiner end
struct EGC <: CorPol end

struct MRC_PSD <: CorPol 
    fs::Int
end

MRC_PSD() = MRC_PSD(10)

"""
    get_weights(x, ref, comb)

Returns the weights for a section of the PSM by the method comb
"""
function get_weights(comb::PCC, x::AbstractMatrix, ref::Int)
    return cor(x, view(x, :, ref)) |> vec
end

function get_weights(comb::PCC2, x::AbstractMatrix, ref::Int)
    pows = sum(abs2, x, dims=1) .|> sqrt
    cors = @views sum(x .* x[:, ref], dims=1)
    w = cors ./ (pows .* pows[ref]) |> vec
    return w
end

function get_weights(comb::SNR_MAX, x::AbstractMatrix, ref::Int)
    cors = @views sum(x .* x[:, ref], dims=1)
    w = cors ./ cors[ref] |> vec
    return w
end

get_weights(comb::Combiner, x::AbstractMatrix) = get_weights(comb, x, choose_ref(x))

function get_weights(comb::MRC_PSD, x::AbstractMatrix{T}) where T
    w = Vector{T}(undef, size(x, 2))
    for (i, s) in enumerate(eachcol(x))
        w[i] = estimate_snr(s, fs=comb.fs)
    end
    return w
end

"""
    combiner(x, comb)

Gets weights from comb and applies them to x
"""
function combiner(comb::Combiner, x::AbstractMatrix)
    w = get_weights(comb, x)
    return x * w
end

function combiner(comb::CorPol, x::AbstractMatrix)
    w = get_weights(comb, x)
    return polarity_flip(x) * w
end

combiner(comb::EGC, x::AbstractMatrix) = sum(polarity_flip(x), dims=2) |> vec

combiner(x::AbstractMatrix) = combiner(SNR_MAX(), x)