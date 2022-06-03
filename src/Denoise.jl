module Denoise
export denoise!, bandpass

import DSP: Highpass, Lowpass, iirnotch, FIRWindow, digitalfilter, filtfilt, kaiser, blackman
include("Raw.jl")
using ..Raw: RawEEG

"""
    `highpass(data::AbstractArray{<:Real}, sampfreq::Int, low::Real)`

Apply highpass filtering to data.

# Return 
An array having the same size as the input array.
"""
function highpass(data::AbstractArray{<:Real}, sampfreq::Int, low::Real)
    wd = 2.0 * low / sampfreq
    filter = digitalfilter(Highpass(wd), 
             FIRWindow(blackman(trunc(Int, sampfreq/2-1))))
    filtfilt(filter, data)
end

"""
    `lowpass(data::AbstractArray{<:Real}, sampfreq::Int, high::Real)`

Apply lowpass filtering to data.

# Return 
An array having the same size as the input array.
"""
function lowpass(data::AbstractArray{<:Real}, sampfreq::Int, high::Real)
    wd = high / sampfreq
    filter = digitalfilter(Lowpass(wd), 
             FIRWindow(blackman(trunc(Int, sampfreq/high - 1))))
    filtfilt(filter, data)
end

"""
    `linenoise_filter(data::AbstractArray{<:Real}, sampfreq::Int, 
                     noisefreq::Real, bandwidth::Real)`

Apply linenoise filtering to data. 

# Return 
An array having the same size as the input array.
"""
function linenoise_filter(
    data::AbstractArray{<:Real}, sampfreq::Int, noisefreq::Real, bandwidth::Real
)
    filter = iirnotch(noisefreq, bandwidth, fs=sampfreq)
    filtfilt(filter, data)
end

"""
    `bandpass(data::AbstractArray{<:Real}, sampfreq::Int; low::Real=1, high::Real=50)`

Apply bandpass filtering to data. 

# Arguments
- `data::AbstractArray{<:Real}`: an array of size (n_channels, n_time)
- `sampfreq::Int`: sampling frequency
- `low::Real=1`: lowest frequency
- `high::Real=50`: highest frequency

# Return 
An array having the same size as the input array.
"""
function bandpass(data::AbstractArray{<:Real}, sampfreq::Int; 
    low::Real=1, high::Real=50
)
    data = transpose(data)
    data = highpass(data, sampfreq, low)
    data = lowpass(data, sampfreq, high)
    transpose(data)
end

"""
    `denoise!(eeg::RawEEG; low::Real=0.8, linenoise::Real=50.0, 
              notch_bandwidth::Real=5.0)`

Denoise data. 

# Arguments
- `eeg::RawEEG`: raw EEG struct
- `low::Real=1`: lowest frequency
- `high::Real=50`: highest frequency
- `notch_bandwidth::Real=5.0`: bandwidth for notch filtering.

# Return 
nothing.

**NOTE: This function causes side effect by modifying `eeg::RawEEG`.**
"""
function denoise!(
    eeg::RawEEG;
    low::Real=0.8, linenoise::Real=50.0, notch_bandwidth::Real=5.0
)
    data = eeg.data
    high = eeg.sampfreq ÷ 2
    data_bp = bandpass(data, eeg.sampfreq; low=low, high=high)
    data_ln = linenoise_filter(transpose(data_bp), eeg.sampfreq, linenoise, notch_bandwidth)
    eeg.data = transpose(data_ln)
end

end