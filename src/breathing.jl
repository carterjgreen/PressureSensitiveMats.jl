function est_br(x::AbstractVector, fs=10)
    # Returns breathing rate in bpm
    bounds = Int.((1.2 * fs, 30 * fs))
    s = autocor(x, bounds[1]:bounds[2])
    delay = argmax(s) + bounds[1]
    return 60 * fs / delay
end

function est_br_fft(x::AbstractVector, fs=10)
    # Returns breathing rate in bpm
    n = length(x)
    s = abs2.(rfft(x))
    f = rfftfreq(n, fs)
    return f[argmax(s)] * 60
end

function est_br_fft2(x::AbstractVector, fs=10)
    # using periodogram, not accurate
    s = periodogram(x, fs=fs)
    return freq(s)[argmax(power(s))] * 60
end