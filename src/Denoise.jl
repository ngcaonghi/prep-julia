module Denoise
export denoise!, bandpass

import DSP: Highpass, Lowpass, iirnotch, FIRWindow, digitalfilter, filtfilt, kaiser, blackman
include("Raw.jl")
using ..Raw: RawEEG

function highpass(data::AbstractArray{<:Real}, sampfreq::Int, low::Real)
    wd = 2.0 * low / sampfreq
    filter = digitalfilter(Highpass(wd), 
             FIRWindow(blackman(trunc(Int, sampfreq/2-1))))
    filtfilt(filter, data)
end

function lowpass(data::AbstractArray{<:Real}, sampfreq::Int, high::Real)
    wd = high / sampfreq
    filter = digitalfilter(Lowpass(wd), 
             FIRWindow(blackman(trunc(Int, sampfreq/high - 1))))
    filtfilt(filter, data)
end

function linenoise_filter(
    data::AbstractArray{<:Real}, sampfreq::Int, noisefreq::Real, bandwidth::Real
)
    filter = iirnotch(noisefreq, bandwidth, fs=sampfreq)
    filtfilt(filter, data)
end

function bandpass(data::AbstractArray{<:Real}, sampfreq::Int; 
    low::Real=1, high::Real=50
)
    data = transpose(data)
    data = highpass(data, sampfreq, low)
    data = lowpass(data, sampfreq, high)
    transpose(data)
end

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