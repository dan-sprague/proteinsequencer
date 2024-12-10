# README 

Instructions for Mac/Linux. 

0. Ensure that some version of conda is installed on your system. This can be done by installing an Anaconda distribution

https://www.anaconda.com/download/success

and BLAST, available via Homebrew for Mac

`brew install blast`

or on Ubuntu

`sudo apt install ncbi-blast+`



1. Install PASS w/ conda 

`conda install bioconda::pass`

2. Install julia (ensure it is on your PATH)

`curl -fsSL https://install.julialang.org | sh`

3. change directory to the GlyphicSequencer folder/git 

`cd GlyphicSequencer/`

4. Run the build script, change the number of cpus to match your computer.

`julia --threads 8 build.jl`


5. Run the simulation(s)

A single simulation


```bash
julia runsim.jl -n $n -a $a -c 0.9 -o 0.94 -l 0.85 -v 0.975 --length 1000 --digest < $PROTEOME_FASTA > $r $READS < $EXP_NAME >
```

For a series of simulations:

```bash
for a in 0.8 0.9 0.95;do for r in "1,1" "1,100" "1,1000"; for n in 30; do for it in {1..3}; do julia runsim.jl -n $n -a $a -c 0.9 -o 0.94 -l 0.85 -v 0.975 --length 1000 --digest < $PROTEOME_FASTA > $r 100000 < $EXP_NAME >;done;done;done
```

    1. n is the number of cycles
    2. a is the base caller accuracy
    3. r is the ratio, can be any number of proteins
    4. it is the number of simulations to run for a set of conditions
    5. positional args are the proteome fasta file and the number of reads


This will create a folder for each set of conditions. The folder name will contain all of the identifying information for the simulation (parameter settings). Inside will be a series of PASS output files.


6. Compile contigs for each set of simulations

```bash
for file in n*;do cat $file/*contigs | sed '/^>/N;/\n[A-Z]\{50,\}$/!d' > $file/sim_contigs.fasta;done
```

7. Run BLAST analysis

```bash
julia scripts/blastscoring.jl $EXP_NAME/
```









