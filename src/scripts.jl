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

    seqs = sample(parse_proteome(args["proteome"];
    minl = args["length"]),n_prots;replace=false)

    for replicate ∈ 1:args["replicates"]
    
        PROTEASES = [r"(?=[R])",r"(?=[K])",r"(?=[E])",r"(?=[M])"]
        peptides = Peptide[]
        if args["digest"]
            for i in 1:n_prots
                seq = seqs[i]
                for protease in PROTEASES
                    nreads = prot_reads[i] ÷ size(PROTEASES,1)
                    for j in 1:nreads
                        frags = digest(seq.sequence,protease)
                        frags = frags[.!(isempty.([i.sequence for i in frags]))]
                        push!(peptides,sample(frags))
                    end
                end
            end
        else
            for i in 1:n_prots
                nreads = prot_reads[i]
                seq = seqs[i]
                for j in 1:nreads
                    push!(peptides,fragment(seq.sequence,λ))
                end
            end
        end

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


        redirect_stdio(stdout=joinpath(OUTDIR,"stdout.txt"), stderr=joinpath(args["outdir"],dir,"stderr.txt")) do
            run(`pass -f $(joinpath(OUTDIR,"reads_$replicate.fasta")) -m 4 -o 1 -r 0.51 -t 0 -w 1 -q 0 -y 0`)
        end
    end

end


function blast_analysis(simsdir)


    p = r".contigs"
    sims = [joinpath(simsdir,i) for i in readdir(simsdir) if !isnothing(match(p,i))]
        
    open("$(simsdir)/sim_contigs.fasta","a") do writer
        for sim in sims
            FASTAReader(open(sim)) do reader
                for record in reader
                    if length(sequence(record)) > 50
                        write(writer,">$(description(record))\n")
                        write(writer,"$(sequence(record))\n")
                    end
                end
            end
        end
    end

    run(`blastp -db /Users/dansprague/Documents/blastdb/swissprot -num_threads 8 -query $(joinpath(simsdir,"sim_contigs.fasta")) -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen staxids" -out $(joinpath(simsdir,"blast_results.tsv"))`)


end


function report()

    p = r"^n\d+"

    
    out = []
    for d in sims
        @show d
        df = CSV.read(joinpath(simsdir,d,"blast_results.tsv"),DataFrame;header=false)
        if isempty(df)
            FASTAReader(open(joinpath(simsdir,d,"sequences.fasta"))) do reader
                
                for record in reader
        
                    name = identifier(record)
                    ratio = (split(name,'|'))[end]
                    name = String.(split(name,'|')[1:2])
                    name = join(name,'|') * '|'
                    push!(out,vcat([i[2:end] for i in split(d,'_')],name,ratio,0,0,0,0,0,0))
                end
        
            end
    
        else 
            seqs = Dict()
            FASTAReader(open(joinpath(simsdir,d,"sequences.fasta"))) do reader
                for record in reader
        
                    name = identifier(record)
                    ratio = (split(name,'|'))[end]
                    name = String.(split(name,'|')[1:2])
                    name = join(name,'|') * '|'
                    seqs[name] = (sequence(record),ratio)
                end
        
            end
    
            
    
            for (name,data) in seqs
                seq,ratio = data
                pos = Set()
                ids = df.Column2
                ids = [i[1] * '|' for i in split.(ids,'.')]
                idx = ids .== name .&& df[:,11] .<= 1e-1
    
                counts = df[idx,:]
    
                taxid = if !isempty(counts) 
                    mean((combine(groupby(counts, 1), [12, 15] => (x, y)->
                y[x .== maximum(x)] .== "6239"))[:,2])
                else
                    0
                end
    
                
                if isempty(counts)
                    push!(out,vcat([i[2:end] for i in split(d,'_')],name,ratio,0,0,0,0,0,0))
                else
                    for row in eachrow(counts)
                        map(x -> push!(pos,x),row[9]:row[10])
                    end
                
                    push!(out,vcat([i[2:end] for i in split(d,'_')],name,round(parse(Float32,ratio);digits=6),maximum(counts[:,13]),maximum(counts[:,14]),maximum(counts[:,13]) ./ maximum(counts[:,14]),length(pos) / counts[1,14],first(counts[counts[:,13] .== maximum(counts[:,13]),3]),taxid))
                end
    
            end
    
            
        end
     
        
    end
    
    out = permutedims(DataFrame(reduce(hcat,out),:auto))
    rename!(out,names(out) .=> [:cycles, :nProts,:acc,:fragL,:ratio,:click,:oligo,:ligate,:cleave,:name,:proportionOfReads,:maxContigLength,:referenceLength,:fractionMaxContigCoverage,:totalCoverage,:maxContigPercIden,:maxBlastTaxIdFracCorrect])
    CSV.write("$(joinpath(ARGS[1],simulation_results.csv))",out)
end
