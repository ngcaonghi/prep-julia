include("../src/Raw.jl")
include("../src/Denoise.jl")
include("../src/Detect.jl")
using Base, Glob, ..Raw, ..Denoise, ..Detect

function process_file(f::AbstractString)
    eeg = loadfile(f)
    println("Starting PREP preprocessing for " * eeg.name * "\n")
    println("Denoising...\n")
    denoise!(eeg)
    println("Looking for bad channels...\n")
    findbad!(eeg)
    println("Bad channels report for " * eeg.name * "\n")
    for (k, v) in eeg.bad
        print("\tBad by " * k *": \t")
        bad = eeg.channels[eeg.usablechans][v]
        if length(bad) == 0
            bad = nothing
        end
        print(bad)
        println()
    end
    println("\nEnd PREP preprocessing for " * eeg.name * "\n===================\n")
end

function process_file()
    println("Input folder directory to .csv file:")
    f = readline()
    process_file(f)
end

function process_folder()
    println("Input folder directory to .csv files:")
    folder = readline()
    pattern = "*.csv"
    files = glob(pattern, folder)
    for f in files
        process_file(f)
    end
end

function main()
    while true
        println("If single file, type \"1\". If whole folder, type \"2\".")
        option = readline()
        if option == "1" 
            process_file()
            break
        elseif option == "2"
            process_folder()
            break
        else
            println("Please choose either 1 (for processing a single file) or 2 (for processing a whole folder).\n")
        end
    end
end

main()