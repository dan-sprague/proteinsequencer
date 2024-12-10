using GlyphicSequencer

"""
    run_simulation()

This scripts runs a simulation given a set of parameters, creates a new directory, and saves the output inside that directory. The name of the directory contains the identifying parameters of the simulation.

CLI arguments are found in GlyphicSequencer.jl

"""
function run_simulation(args,OUTDIR)
    prot_reads = Int.(ceil.(args["reads"]))

    isdir(args["outdir"]) ? nothing : mkdir(args["outdir"])

    println("###### LOG INFO #######")
    @show args["ncycles"]
    @show args["reads"]


    @show prot_reads
    @show args["replicates"]

    @show args["acc"]

    sequencer = Sequencer(args["click"],args["oligo"],args["ligate"],args["cleave"],BaseCaller(args["acc"],AA))
    
    @show sequencer

    seqs = parse_proteome(args["peptides"];minl = 9)
    peptides = [Peptide(seq.sequence,'E','C','O',"",1) for seq in seqs]


    @assert all(i.sequence != "" for i in peptides)
    
    @inbounds for i ∈ 1:NCYCLES
        @. peptides = sequencer(peptides)
    end
    
    @show length(peptides)

    open(joinpath(OUTDIR,"sequences.fasta"),"a") do file
        i = 1
        for seq in seqs
            write(file,">$(seq.id)\n")
            write(file,"$(seq.sequence)]\n")
            i += 1
        end
    end


    open(joinpath(OUTDIR,"reads.fasta"),"w") do file
        i = 0
        for read in peptides
            write(file,">$i\n")
            write(file,"$(read.aas)\n")
            i += 1
        end
    end

end


args = parse_commandline()
dir = "n$(args["ncycles"])_a$(args["acc"])_c$(args["click"])_o$(args["oligo"])_l$(args["ligate"])_v$(args["cleave"])"
isdir(joinpath(args["outdir"],dir)) ? nothing : mkdir(joinpath(args["outdir"],dir)) 

OUTDIR = joinpath(args["outdir"],dir)
run_simulation(args,OUTDIR)

