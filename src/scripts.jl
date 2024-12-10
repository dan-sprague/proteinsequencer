"""
    run_simulation()

This scripts runs a simulation given a set of parameters, creates a new directory, and saves the output inside that directory. The name of the directory contains the identifying parameters of the simulation.

CLI arguments are found in GlyphicSequencer.jl

"""
function run_simulation(args,OUTDIR)
    ratio = parse.(Int,split(args["ratio"],","))
    r_norm = ratio ./ sum(ratio)
    prot_reads = Int.(ceil.(args["reads"] .* r_norm))
    n_prots = size(ratio,1)
    λ = args["lambda"]

    isdir(args["outdir"]) ? nothing : mkdir(args["outdir"])

    println("###### LOG INFO #######")
    @show args["ncycles"]
    @show args["reads"]

    @show n_prots
    @show ratio
    @show r_norm
    @show prot_reads
    @show λ
    @show args["replicates"]

    @show args["acc"]

    sequencer = Sequencer(args["click"],args["oligo"],args["ligate"],args["cleave"],BaseCaller(args["acc"],AA))
    
    @show sequencer

    seqs = parse_proteome(args["proteome"];minl = 9)
    peptides = [Peptide(seq,'E','C','O',"",1) for seq in seqs]

    for replicate ∈ 1:args["replicates"]
    
        peptides = Peptide[]

        @assert all(i.sequence != "" for i in peptides)


        μ_length = Vector{Float32}(undef,args["ncycles"])
        σ_length = Vector{Float32}(undef,args["ncycles"])

        μ_pos = Vector{Float32}(undef,args["ncycles"])
        σ_pos = Vector{Float32}(undef,args["ncycles"])

        μ_deletions = Vector{Float32}(undef,args["ncycles"])
        σ_deletions = Vector{Float32}(undef,args["ncycles"])

        dels  = Vector{Int}(undef,length(peptides))
        X = (μ_length,σ_length,μ_pos,σ_pos,μ_deletions,σ_deletions,dels)
        
        simulate!(sequencer,peptides,X;NCYCLES = args["ncycles"])
        

        open(joinpath(OUTDIR,"sequences.fasta"),"a") do file
            i = 1
            for seq in seqs
                write(file,">$(seq.id)|$(r_norm[i])\n")
                write(file,"$(seq.sequence)]\n")
                i += 1
            end
        end


        open(joinpath(OUTDIR,"reads_$(replicate).fasta"),"w") do file
            i = 0
            for read in peptides
                write(file,">$i\n")
                write(file,"$(read.aas)\n")
                i += 1
            end
        end

    end

end