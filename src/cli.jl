function parse_commandline()
    s = ArgParseSettings()

    @add_arg_table s begin
        "--ncycles","-n"
            help = "Cycles of sequencing"
            arg_type = Int
            default = 15
        
        "--click","-c"
            help = "Click probability"
            arg_type = Float32
            default = .9
        "--oligo","-o"
            help = "Oligo probability"
            arg_type = Float32
            default = .9
        "--ligate","-l"
            help = "Ligate probability"
            arg_type = Float32
            default = .9
        "--cleave","-v"
            help = "Cleave probability"
            arg_type = Float32
            default = .9

        "--loss","-d"
            help = "Probability of sequence loss"
            arg_type = Float32
            default = .9
        "--acc","-a"
            help = "Average basecalling accuracy"
            arg_type = Float32
            default = 0.8
            
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