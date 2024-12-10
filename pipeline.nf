#!/usr/bin/env nextflow

params.proteome = "$baseDir/proteome.fasta"
params.reads = 100000
params.lambda = 15
params.cycles = 15
params.ratio = "1,1"
params.outdir = "sim"
params.click = 0.9
params.oligo = 0.9
params.ligate = 0.9
params.cleave = 0.9
params.acc = 0.9
params.length = 1000

/*
 * Split a fasta file into multiple files
 */
process simulation {

    input:
    val click
    val oligo
    val ligate
    val cleave
    val acc
    val length
    path proteome
    val ratio
    val reads
    val lambda
    val cycles
    val outdir



    output:
    path outdir

    """
    julia $baseDir/scripts/runsim.jl -n $cycles -c $click -o $oligo -l $ligate -v $cleave -a $acc --length $length --lambda $lambda $proteome $ratio $reads $outdir
    """
}


process prepare_for_blast {

    input:
    path file

    output:
    path file


    """
    cat $file/*contigs > $file/sim_contigs.fasta
    julia $baseDir/parse_contigs.jl $file/sim_contigs.fasta
    """

}

/*
 * Define the workflow
 */
workflow {

    simulation(params.click,params.oligo,params.ligate,params.cleave,params.acc,params.length,params.proteome,params.ratio,params.reads,params.lambda,params.cycles,params.outdir)

    dirs = Channel
            .fromPath("${params.outdir}/n*")
            .view()


}