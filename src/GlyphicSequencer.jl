module GlyphicSequencer

using FASTX

using StatsBase: sample,mean,std
using ArgParse
using DataStructures 
using Glob 


include("preprocess.jl")
include("peptide.jl")
include("basecaller.jl")
include("sequencer.jl")
include("cli.jl")


global const AA = collect("ACDEFGHIKLMNPQRSTVWY")


export AA 
export parse_commandline
export Sequencer 

end # module GlyphicSequencer