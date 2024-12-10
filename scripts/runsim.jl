using GlyphicSequencer

args = parse_commandline()
ratio = parse.(Int,split(args["ratio"],","))

ratio_var = join(ratio,'.')
dir = "n$(args["ncycles"])_p$(n_prots)_a$(args["acc"])_f$(args["lambda"])_r$(ratio_var)_c$(args["click"])_o$(args["oligo"])_l$(args["ligate"])_v$(args["cleave"])"
isdir(joinpath(args["outdir"],dir)) ? nothing : mkdir(joinpath(args["outdir"],dir)) 

OUTDIR = joinpath(args["outdir"],dir)
run_simulation(args,OUTDIR)
blast_analysis(outdir)

