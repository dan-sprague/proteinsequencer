module GlyphicSequencer

using FASTX
using Plots 
using StatsBase: sample,mean,std
using Distributions: Poisson
using ArgParse
using DataStructures 
using Glob 
using DataFrames
using CSV 


include("preprocess.jl")
include("peptide.jl")
include("basecaller.jl")
include("sequencer.jl")
include("simulation.jl")
include("scripts.jl")



global const AA = collect("ACDEFGHIKLMNPQRSTVWY")

function parse_commandline()
    s = ArgParseSettings()

    @add_arg_table s begin
        "--ncycles","-n"
            help = "Path to proteome fasta file"
            arg_type = Int
            default = 15
        
        "--click","-c"
            help = "Path to proteome fasta file"
            arg_type = Float32
            default = .9
        "--oligo","-o"
            help = "Path to proteome fasta file"
            arg_type = Float32
            default = .9
        "--ligate","-l"
            help = "Path to proteome fasta file"
            arg_type = Float32
            default = .9
        "--cleave","-v"
            help = "Path to proteome fasta file"
            arg_type = Float32
            default = .9
        "--acc","-a"
            help = "Average basecalling accuracy"
            arg_type = Float32
            default = 0.8


        "--replicates"
            help = "Number of simulations"
            arg_type = Int
            default = 5

        "peptides"
            help = "Path to fasta file"
            required = true
        
        "reads"
            help = "Number of reads"
            arg_type = Int
            required = true

        "outdir"
            help = "Output directory name"
            arg_type = String
            required = true 
    end

    return parse_args(s)
end

export parse_commandline
export run_simulation
export blast_analysis
export parse_contigs



end # module GlyphicSequencer