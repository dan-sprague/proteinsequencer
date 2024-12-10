using Distributions
using GlyphicSequencer
conditions = round.(rand(Beta(5,1),4,300);digits=4)

args = parse_commandline()
ratio = 1
r_norm = ratio ./ sum(ratio)
prot_reads = Int.(ceil.(args["reads"] .* r_norm))
n_prots = 1
λ = 15

println("###### LOG INFO #######")
@show args["ncycles"]
@show args["reads"]

@show n_prots
@show ratio
@show r_norm
@show prot_reads
@show λ

seqs = sample(parse_proteome(args["proteome"];
            minl = args["length"]),n_prots;replace=false)
isdir(args["outdir"]) ? nothing : mkdir(args["outdir"]) 
cd(args["outdir"])


sims = [Sequencer(x[1:4]...,BaseCaller(rand(Beta(5,1)),AA)) for x in eachcol(conditions)]


peptides = Peptide[]
for i in 1:n_prots
    nreads = prot_reads[i]
    seq = seqs[i]
    for j in 1:nreads
        push!(peptides,fragment(seq.sequence,λ))
    end
end


μ_length = Vector{Float32}(undef,args["ncycles"])
σ_length = Vector{Float32}(undef,args["ncycles"])

μ_pos = Vector{Float32}(undef,args["ncycles"])
σ_pos = Vector{Float32}(undef,args["ncycles"])

μ_deletions = Vector{Float32}(undef,args["ncycles"])
σ_deletions = Vector{Float32}(undef,args["ncycles"])

dels  = Vector{Int}(undef,length(peptides))
X = (μ_length,σ_length,μ_pos,σ_pos,μ_deletions,σ_deletions,dels)


Threads.@threads for sim in sims

    simulate!(sim,peptides,X;NCYCLES = args["ncycles"])
    ratio_var = join(ratio,'.')
    dir = "n$(args["ncycles"])_p$(n_prots)_a$(sim.basecall.acc)_f$(args["lambda"])_r$(ratio_var)_c$(sim.c_rate)_o$(sim.o_rate)_l$(sim.l_rate)_v$(sim.clv_rate)"
    isdir(dir) ? nothing : mkdir(dir) 

    open(joinpath(dir,"sequences.fasta"),"a") do file
        i = 1
        for seq in seqs
            write(file,">$(seq.id)|$(r_norm[i])\n")
            write(file,"$(seq.sequence)]\n")
            i += 1
        end
    end

    open(joinpath(dir,"reads.fasta"),"w") do file
        i = 0
        for read in peptides
            write(file,">$i\n")
            write(file,"$(read.aas)\n")
            i += 1
        end
    end

    println("\n\n\n###### PASS #######\n\n\n")
    redirect_stdio(stdout=joinpath(dir,"stdout.txt"), stderr=joinpath(dir,"stderr.txt")) do
        run(`pass -f $(joinpath(dir,"reads.fasta")) -m 4 -o 1 -r 0.51 -t 0 -w 1 -q 0 -y 0`)
    end


end