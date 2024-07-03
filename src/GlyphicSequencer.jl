module GlyphicSequencer

using FASTX
using Plots 
using StatsBase: sample,mean,std
using Distributions: Poisson
using ArgParse
using DataStructures 
using Glob 
using DataFrames
include("preprocess.jl")
include("peptide.jl")
include("basecaller.jl")
include("sequencer.jl")
include("simulation.jl")
include("passanalysis.jl")

global const AA = collect("ACDEFGHIKLMNPQRSTVWY")

function parse_commandline()
    s = ArgParseSettings()

    @add_arg_table s begin
        "--ncycles","-n"
            help = "Path to proteome fasta file"
            arg_type = Int
            default = 15
        "--ratio", "-r"
            help = "Ratio of peptides in mix, for example \"[1,100,1000]\". This also defines the number of peptides that will be sampled."
            default = "1,1"

        "--acc","-a"
            help = "Average basecalling accuracy"
            arg_type = Float32
            default = 0.8
        "--length","-l"
            help = "Protein length"
            default = 1000
        "proteome"
            help = "Path to fasta file"
            required = true
        "λ"
            help = "Average fragment size"
            arg_type = Int
            required = true            
        "reads"
            help = "Number of reads"
            arg_type = Int
            required = true
    end

    return parse_args(s)
end

export AA 
export Sequencer
export BaseCaller
export Peptide
export simulate!
export plot
export digest
export fragment
export parse_proteome
export parse_commandline
export sample
export pass_report 
end # module GlyphicSequencer