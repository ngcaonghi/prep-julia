include("Raw.jl")
using Raw

mutable struct PrepEEG
    data::RawEEG
    highpass::Real
    lowpass::Real
    linenoise::Array{Real} 
end