"""
	Peptide(sequence::String,state::Char,oligo::Char,ligate::Char,aas::String,pos::Int64)

A data type that contains all the necessary features to track a peptide through a sequencing run.
"""
struct Peptide
	sequence::String
	state::Char
	oligo::Char
	ligate::Char
	aas::String
	pos::Int64
end

"""
    digest(peptide,aa;term = 'C')

Simulates the digestion of a `peptide`` given a certain enzyme defined by `aa` and cleavage location `term`. Returns a vector of Peptide objects ready to be loaded on to the in silico sequencer.
"""
function digest(peptide,aa;term='C')
    if term == 'C'
        j = vcat(0,findall(x -> x == aa,peptide),length(peptide))
        [peptide[j[i]+1:j[i+1]] for i in 1:size(j,1) - 1]
    else
        error("NotImplementedError")
    end
end

"""
    fragment(peptide,λ)

This function generates random fragments of `peptide`. The length of these fragments follow a poisson distribution with mean (and definitionally variance) =  λ.

"""
function fragment(peptide,λ)

	frag_length = rand(Poisson(λ))
	pos = rand(1:length(peptide)-frag_length)
    
	Peptide(peptide[pos:pos+frag_length],'E','C','O',"",1)
end