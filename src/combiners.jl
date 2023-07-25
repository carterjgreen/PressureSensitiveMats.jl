abstract type Combiner end
abstract type CorPol <: Combiner end

struct PCC <: Combiner end
struct PCC2 <: Combiner end
struct SNR_MAX <: Combiner end
struct EGC <: CorPol end

struct MRC_PSD <: CorPol
    fs::Real
end

MRC_PSD() = MRC_PSD(10)

"""
    get_weights(comb, x, ref)

Returns the weights for a section of the PSM by the method comb
"""
function get_weights(comb::PCC, x::AbstractMatrix{<:Number}, ref::Integer)
    return vec(cor(x, view(x, :, ref)))
end

function get_weights(comb::PCC2, x::AbstractMatrix{<:Number}, ref::Integer)
    pows = sqrt.(sum(abs2, x; dims = 1))
    cors = @views sum(x .* x[:, ref], dims = 1)
    w = vec(cors ./ (pows .* maximum(pows)))
    return w
end

function get_weights(comb::SNR_MAX, x::AbstractMatrix{<:Number}, ref::Integer)
    cors = @views sum(x .* x[:, ref], dims = 1)
    w = vec(cors ./ maximum(cors))
    return w
end

function get_weights(comb::Combiner, x::AbstractMatrix{<:Number})
    return get_weights(comb, x, choose_ref(x))
end

function get_weights(comb::MRC_PSD, x::AbstractMatrix{T}) where {T <: Number}
    w = Vector{T}(undef, size(x, 2))
    for (i, s) in enumerate(eachcol(x))
        w[i] = estimate_snr(s; fs = comb.fs)
    end
    return w
end

"""
    combiner(comb, x, ref)

Gets weights from combining method comb with a chosen reference sensor and applies it to x.    
"""
function combiner(comb::Combiner, x::AbstractMatrix{<:Number}, ref::Integer)
    w = get_weights(comb, x, ref)
    return x * w
end
"""
    combiner(comb, x)

Gets weights from combining method comb and applies them to x
"""
function combiner(comb::Combiner, x::AbstractMatrix{<:Number})
    w = get_weights(comb, x)
    return x * w
end

function combiner(comb::CorPol, x::AbstractMatrix{<:Number})
    w = get_weights(comb, x)
    return polarity_flip(x) * w
end

combiner(comb::EGC, x::AbstractMatrix{<:Number}) = vec(sum(polarity_flip(x); dims = 2))

combiner(x::AbstractMatrix{<:Number}) = combiner(SNR_MAX(), x)

"""
    snrmax(x)

Perform SNR-MAX combining on x. Equivalent to combiner(SNR_MAX(), x)
"""
snrmax(x::AbstractMatrix{<:Number}) = combiner(SNR_MAX(), x)

"""
    mrc(x, fs)

Perform MRC-PSD combining on x. Currently x must be less than 512 samples.
Equivalent to combiner(MRC_PSD(fs), x)
"""
mrc(x::AbstractMatrix{<:Number}, fs::Number) = combiner(MRC_PSD(fs), x)

"""
    egc(x)

Perform equal gain combining on x. Equivalent to combiner(EGC(), x)
"""
egc(x::AbstractMatrix{<:Number}) = combiner(EGC(), x)

"""
    pcc(x)

Perform Pearson Correlation Coefficient combining on x.
Equivalent to combiner(PCC2(), x)
"""
pcc(x::AbstractMatrix{<:Number}) = combiner(PCC2(), x)

"""
    selection(x)

Perform selection combining on x.
Equivalent to extract_ref(x)
"""
selection(x::AbstractMatrix{<:Number}) = extract_ref(x)
