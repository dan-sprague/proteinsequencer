using GlyphicSequencer

args = parse_commandline()
dir = "n$(args["ncycles"])_a$(args["acc"])_c$(args["click"])_o$(args["oligo"])_l$(args["ligate"])_v$(args["cleave"])"
isdir(joinpath(args["outdir"],dir)) ? nothing : mkdir(joinpath(args["outdir"],dir)) 

OUTDIR = joinpath(args["outdir"],dir)
run_simulation(args,OUTDIR)

