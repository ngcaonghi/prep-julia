include("../src/Raw.jl")
include("../src/Denoise.jl")
include("../src/Detect.jl")
using Base, Glob, ..Raw, ..Denoise, ..Detect

"""
    process_file(f::AbstractString)

Process file from the input directory `f::AbstractString`.

**NOTE: Only denoising and bad channel detection have been implemented.**
"""
function process_file(f::AbstractString)
    eeg = loadfile(f)
    println("Starting PREP preprocessing for " * eeg.name * "\n")
    println("Denoising...\n")
    denoise!(eeg)
    println("Looking for bad channels...\n")
    findbad!(eeg)
    println("Bad channels report for " * eeg.name * "\n")
    for (k, v) in eeg.bad
        print("\tBad by " * k *": ")
        bad = eeg.channels[eeg.usablechans][v]
        if length(bad) == 0
            bad = nothing
        end
        print(bad)
        println()
    end
    println("\nEnd PREP preprocessing for " * eeg.name * "\n====================================\n")
end

"""
    process_file()

Ask user to provide a file name leading to the file that needs preprocessing.
"""
function process_file()
    println("Input folder directory to .csv file:")
    f = readline()
    process_file(f)
end

"""
    process_file()

Ask user to provide a folder name leading to the files that need preprocessing.
"""
function process_folder()
    println("Input folder directory to .csv files:")
    folder = readline()
    pattern = "*.csv"
    files = glob(pattern, folder)
    for f in files
        process_file(f)
    end
end

"""
    main()

Driver method.
"""
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