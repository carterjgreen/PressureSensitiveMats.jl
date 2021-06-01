abstract type Combiner end
abstract type CorPol <: Combiner end

struct PCC <: Combiner end
struct SNR_MAX <: Combiner end
struct MRC_PSD <: CorPol end
struct EGC <: CorPol end

"""
    get_weights(x, ref, comb)

Returns the weights for a section of the PSM by the method comb
"""
function get_weights(x::AbstractMatrix, ref::Int, comb::PCC; mean=true)
    if mean
        w = w = cor(x, x[:, ref]) |> vec
    else
        pows = sum(abs2, x, dims=1) .|> sqrt
        cors = @views sum(x .* x[:, ref], dims=1)
        w = cors ./ (pows .* pows[ref]) |> vec
    end
    return w
end

function get_weights(x::AbstractMatrix, comb::PCC; mean=true)
    return get_weights(x, choose_ref(x), comb, mean=mean)
end

function get_weights(x::AbstractMatrix, ref::Int, comb::SNR_MAX)
    # Also called MVDR in Hilda and Sareh's work
    cors = @views sum(x .* x[:, ref], dims=1)
    w = cors ./ cors[ref] |> vec
    return w
end

get_weights(x::AbstractMatrix, comb::SNR_MAX) = get_weights(x, choose_ref(x), comb)

function estimate_snr(x::AbstractVector; fs=10)
    # Estimate SNR for a 30-102.4s segment
    pow = power(periodogram(x, nfft=1024, window=hanning))
    ps = argmax(pow)
    n = mean(pow[Not(ps)]) * 0.73
    w = (pow[ps] * fs / 1024 - n) / n
    return w
end

function get_weights(x::AbstractMatrix{T}, comb::MRC_PSD; fs=10) where T
    w = Vector{T}(undef, size(x, 2))
    for (i, s) in enumerate(eachcol(x))
        w[i] = estimate_snr(s, fs)
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
function combiner(x::AbstractMatrix, comb::Combiner, kwargs...)
    w = get_weights(x, comb, kwargs...)
    return x * w
end

function combiner(x::AbstractMatrix, comb::CorPol)
    w = get_weights(x, comb)
    return polarity_flip(x) * w
end

combiner(x::AbstractMatrix) = combiner(x, SNR_MAX())