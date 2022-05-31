module Rereference
include("Raw.jl")
include("Detect.jl")
include("Denoise.jl")
include("ChannelPositions.jl")
using .Raw, .Detect, .Denoise, DSP, Base, Statistics
import .ChannelPositions:positions as eegpos

function init_reref(data::AbstractArray{<:Real}; dims=1)
    data .- median(skipmissing(data); dims=dims)
end


"""
    EEG_temp = EEG - robust estimation of mean (median by default)
    badChannels = []
    iterations = 0
    while true
        detect bad channels using findNoisyChannel
        add newly detected bad channels to array
        break from loop if bad channels did not change or iteration criteria
            had been met
        newMean = mean of EEG with all current bad channels interpolated
        EEGTemp = EEG - newMean
        iterations += 1
    end
    reference = mean of EEG with current bad signals interpolated
"""
function phaseone(eeg::RawEEG, epochs::Int; kwargs...)
    findbad(eeg; kwargs...)
    usablemask = vec(eeg.usablechans .& .!eeg.bad["flat"] .& .!eeg.bad["SNR"])
    initref = init_reref(eeg.data_out[usablemask])
    eeg_temp = freshstart(initref, eeg)
    i = 0
    prevbad = deepcopy(eeg.bad)
    cumbad = cumulativebad(prevbad, size(eeg.channels, 1))
    id_range = collect(1:size(eeg.channels, 1))
    while i < epochs
        findbad(eeg_temp, referencemode=true, filter=false)
        nowcumbad = cumulativebad(eeg_temp.bad, size(eeg_temp.channels, 1))
        if (cumbad == nowcumbad)|(i > 1 & sum(nowcumbad)==sumverybad(eeg_temp))
            break
        end
        from = id_range[.!nowcumbad]
        to = id_range[nowcumbad]
        EEG_interp = interpolate(eeg, from, to)
    end
end

function freshstart(
    data::Array{<:Real},
    template::RawEEG
)
    data_out = deepcopy(data)
    bad = Dict()
    bad["nan"] = template.bad["nan"]
    bad["SNR"] = template.bad["SNR"]
    bad["all"] = template.bad["all"]
    return RawEEG(
        data,
        data_out,
        template.channels,
        bad,
        template.usablechans,
        template.sampfreq,
        template.name,
        template.dir
    )
end

function cumulativebad(eegbad::Dict{AbstractString, AbstractArray}, numchans)
    cumbad = zeros(Bool, numchans)
    for bad in eegbad
        cumbad = cumbad .| bad[2]
    end
    cumbad
end

function sumverybad(eeg::RawEEG)
    sum(eeg.bad["nan"]) + sum(eeg.bad["SNR"]) + sum(eeg.bad["all"])
end

function getpos()
end

function interpolate(eeg::RawEEG, from::AbstractArray, to::AbstractArray)
    gfrom_to = sphere_spline_g(from, to)
    gto_from = sphere_spline_g(to, from)
    avg = mean(eeg.data_filtered, dims=1)
    data_now = eeg.data_filtered .- avg
    padding = ones(Float64, 1, size(from)[1])
    C = pinv(vcat(gfrom_to, padding))
    I = gto_from * C[:, 1:size(C, 1)-1]
    I * data_now .+ avg
end

function sphere_spline_g(from, to; stiffness=4, lterms=7)
    nfrom = size(from, 1)
    nto = size(to, 1)
    dsquareds = []
    for i = 1:3
        dist1 = reshape(repeat(to[:, i], outer=nfrom), (nto, nfrom))
        dist2 = transpose(reshape(repeat(from[:, i], outer=nto), (nfrom, nto)))
        push!(dsquareds, (dist1 .- dist2).^2)
    end
    distances = sqrt.(sum(dsquareds))
    EI = ones(Float64, nto, nfrom) .- distances
    factors = [0.0]
    for n = 1:lterms
        f = (2 * n + 1) / (n^stiffness * (n + 1)^stiffness * 4 * π)
        push!(f, n)
    end
    legandre(EI, factors)
end

function legandre(EI::AbstractArray, factors::AbstractArray)
    factors = vec(factors)
    if size(factors, 1) == 1
        f1 = factors[1]
        f2 = 0
    elseif size(factors, 1) == 2
        f1 = factors[1]
        f2 = factors[2]
    else
        n = size(factors, 1)
        f1 = factors[end-1]
        f2 = factors[end]
        for i in 2:size(factors, 1)-1
            tmp = f1
            n = n - 1
            f1 = factors[end-i] .- (f2 .* (n - 1)) / n
            println(f1)
            f2 = tmp .+ (f2.*EI.*(2*n-1)) ./ n
        end
    end 
    f1 .+ f2 .* EI
end

end