module Raw
export RawEEG, loadfile, getchannel
import CSV.read as readcsv
using Base, DataFrames

"""
    RawEEG

# Attributes
- `data::AbstractArray{<:Real}`: data with size (n_channels, n_time)
- `data_filtered::AbstractArray{<:Real}`: bandpass-filtered data
- `channels::AbstractArray{String}`: array of channel names
- `bad::DictDict{AbstractString, AbstractArray}`: dictionary of bad type (String)
=> Bool mask (Array of size (n_channels))
- `usablechans::BitVector`: Bool mask of channels that are not bad by NaN.
- `sampfreq::Int`: sampling frequency
- `name::String`: sample name
- `dir::String`: local file directory
"""
mutable struct RawEEG
    data::AbstractArray{<:Real}
    data_filtered::AbstractArray{<:Real}
    channels::AbstractArray{String}
    bad::Dict{AbstractString, AbstractArray}
    usablechans::BitVector
    sampfreq::Int
    name::String
    dir::String
end
RawEEG(dir; kwargs...) = loadfile(dir; kwargs...)

"""
    loadfile(dir::String;
    channels::Array{String} = ["F3";"FC5";"AF3";"F7";"T7";"P7";"O1";"O2";"P8";"T8";"F8";"AF4";"FC6";"F4"],
    bad::Dict = Dict{AbstractString, AbstractArray}(),
    usablechans::BitVector = vec(channels .!= nothing),
    sampfreq::Int = 256)

Create RawEEG struct from file directory.
"""
function loadfile(
    dir::String;
    channels::Array{String} = [
        "F3";"FC5";"AF3";"F7";"T7";"P7";"O1";"O2";"P8";"T8";"F8";"AF4";"FC6";"F4"
    ],
    bad::Dict = Dict{AbstractString, AbstractArray}(),
    usablechans::BitVector = vec(channels .!= nothing),
    sampfreq::Int = 256
)
    arr_raw = readcsv(dir, DataFrame, header=false, skipto=2)[:, 3:length(channels)+2]
    r = size(arr_raw, 1)
    c = size(arr_raw, 2)
    arr = zeros(Float64, c, r)
    for i in 1:c
        for j in 1:r
            arr[i, j] = arr_raw[j, i]
        end
    end 
    data_filtered = deepcopy(arr)
    name = splitext(basename(dir))[1]
    return RawEEG(
        arr,
        data_filtered,
        channels,
        bad,
        usablechans,
        sampfreq,
        name,
        dir
    )
end

"""
    getchannel(eeg::RawEEG, channel::AbstractString)

Pull out data of a channel by name.
"""
function getchannel(eeg::RawEEG, channel::AbstractString) 
    arr = eeg.data
    idxarr = arr[arr.==channel]
    re_arr = arr[idxarr[end]]
    if isempty(re_arr)
        throw(UndefVarError("No " * channel * " channel in this data."))
    else
        return re_arr
    end
end

end