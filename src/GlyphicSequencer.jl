module GlyphicSequencer

using FASTX

using StatsBase: sample,mean,std
using ArgParse
using Glob 


include("preprocess.jl")
include("peptide.jl")
include("basecaller.jl")
include("sequencer.jl")
include("cli.jl")

global const AA = collect("ACDEFGHIKLMNPQRSTVWY")


export parse_commandline
export AA 
export Sequencer 
export BaseCaller
export parse_proteome
export Peptide

end # module GlyphicSequencer