module Detect
include("Raw.jl")
include("Denoise.jl")
using Base, DSP, Statistics, LinearAlgebra
import ..Raw: RawEEG
import ..Denoise: denoise!, bandpass
export findbad!

function findbad!(
    eeg::RawEEG;
    flatthresh::Real=1e-15,
    devthresh::Real=5.0,
    hfthresh::Real=5.0,
    corrsecs::Real=1.0,
    corrthresh::Real=0.4,
    frac::Real=0.01,
    referencemode=false,
    filter=true
    )
    if referencemode==false
        badnan!(eeg)
    end
    badflat!(eeg; flatthresh=flatthresh)
    baddev!(eeg; devthresh=devthresh)
    if filter==true
        filtered = bandpass(eeg.data_filtered, eeg.sampfreq; low=1, high=50)
        eeg.data_filtered = filtered
    end
    badhf!(eeg; hfthresh=hfthresh)
    badcorr!(eeg; corrsecs= corrsecs, corrthresh=corrthresh, frac=frac)
    if referencemode==false
        badSNR!(eeg)
        badall!(eeg)
    end
end

function badnan!(eeg::RawEEG)
    arr = eeg.data
    sumarr = sum(arr, dims=2) 
    nanmask = sumarr .== NaN
    eeg.bad["nan"] = vec(nanmask)
    eeg.usablechans = vec(.!nanmask)
end

function badflat!(eeg::RawEEG; flatthresh::Real=1e-15)
    arr = eeg.data
    madflat = medabsdev(arr) .< flatthresh
    stdevflat = std(arr, dims=2) .< flatthresh
    flatmask = madflat .| stdevflat 
    eeg.bad["flat"] = vec(flatmask)
end

function baddev!(eeg::RawEEG; devthresh::Real=5.0)
    IQR_TO_SD = 0.7413
    channel_amp = iqr(eeg.data) .* IQR_TO_SD
    channel_amp = vec(channel_amp[eeg.usablechans, :])
    ampsd = (quantile(channel_amp, 0.75) - quantile(channel_amp, 0.25)) * IQR_TO_SD
    ampmed = median(skipmissing(channel_amp))
    ampzscore = (channel_amp .- ampmed) ./ ampsd
    devmask = (ampzscore .== NaN) .| (abs.(ampzscore) .> devthresh)
    eeg.bad["deviation"] = vec(devmask)
end

function badhf!(eeg::RawEEG; hfthresh::Real=5.0)
    MEDABS_TO_SD = 1.4826 
    usable = eeg.data[eeg.usablechans, :]
    filtered = eeg.data_filtered[eeg.usablechans, :]
    noisiness = medabsdev(usable .- filtered) ./ medabsdev(filtered)
    noisemed = median(skipmissing(noisiness))
    noisesd = median(abs.(noisiness .- noisemed)) * MEDABS_TO_SD 
    zscore = (noisiness .- noisemed) ./ noisesd
    hfmask = (zscore .== NaN) .| (zscore .> hfthresh)
    eeg.bad["high_frequency"] = vec(hfmask)
end

function badcorr!(eeg::RawEEG; 
    corrsecs::Real=1.0, corrthresh::Real=0.4, frac::Real=0.01
    )
    winsize = trunc(Int64, corrsecs * eeg.sampfreq)
    winoffsets = collect(1:winsize:size(eeg.data, 1)-winsize)
    wincount = size(winoffsets, 1)
    maxcorr = ones(wincount, size(eeg.usablechans, 1))
    dropout = zeros(Bool, wincount, size(eeg.usablechans, 1))
    for w = 1:wincount
        start = w * winsize
        stop = (w + 1) * winsize
        filtered = eeg.data_filtered[eeg.usablechans, start:stop]
        raw = eeg.data[eeg.usablechans, start:stop]
        usable = deepcopy(eeg.usablechans)
        amp = medabsdev(filtered)
        dropout[w, usable] = amp .== 0
        usable[usable] = amp .> 0
        raw = raw[amp .> 0, :]
        filtered = filtered[amp .> 0, :]
        wincorr = cor(filtered, dims=2)
        abscorr = abs.(wincorr .- diag(diag(wincorr)))
        qt = zeros(size(abscorr, 1))
        for i = 1:size(qt, 1)
            qt[i] = quantile(abscorr[i,:], 0.98)
        end
        maxcorr[w, usable] = qt
        maxcorr[w, dropout[w, :]] .= 0
    end
    threshedcorr = maxcorr .< corrthresh
    fracbad = mean(threshedcorr, dims=1)
    corrmask = fracbad .> frac
    eeg.bad["correlation"] = vec(corrmask)
    fracdrop = mean(dropout, dims=1)
    dropmask = fracdrop .> frac
    eeg.bad["dropout"] = vec(dropmask)
end

function badSNR!(eeg::RawEEG)
    eeg.bad["SNR"] = eeg.bad["correlation"] .& eeg.bad["high_frequency"]
end

function badall!(eeg::RawEEG)
    allmask = zeros(Bool, size(eeg.usablechans, 1))
    for bad in eeg.bad
        if bad[1] == "nan"
            continue
        else
            allmask = allmask .& bad[2]
        end
    end
    eeg.bad["all"] = allmask
end

function iqr(arr::AbstractArray{<:Real})
    qt75 = zeros(Float64, size(arr,1))
    qt25 = zeros(Float64, size(arr,1))
    for i = 1:size(arr,1)
        qt75[i] = quantile(arr[i,:], 0.75)
        qt25[i] = quantile(arr[i,:], 0.25)
    end
    qt75 .- qt25
end

function medabsdev(arr::AbstractArray{<:Real}; dims::Int=2)
    med = median(arr, dims=dims)
    median(abs.(arr .- med), dims=dims)
end

end