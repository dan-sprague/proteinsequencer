"""
    Sequencer(c_rate::Float32,o_rate::Float32,l_rate::Float32,clv_rate::Float32)

An in silico sequencer with tunable properties.

"""
struct Sequencer
    c_rate::Float32
    o_rate::Float32
    l_rate::Float32
    clv_rate::Float32
    basecall::BaseCaller

    function Sequencer(a,b,c,d,e)
        @assert all([a,b,c,d] .<= 1) && all([a,b,c,d] .>= 0)
        new(a,b,c,d,e)
    end
end



"""
    (s::Sequencer)(peptide)

The sequencer `s` is called as a function that takes in the state for a given peptide and applies the sequencing protocol with parameters defined in the in silico sequencer s.

```julia-repl
sequencer = Sequencer(0.85,0.90,0.85,0.85)

sequencer(peptide)

```

"""
function (s::Sequencer)(x::Peptide)
    
    function click(x::Char)
        rand() < s.c_rate && x == 'E' ? 'C' : x
    end

    function oligo(x::Char)
        rand() < s.o_rate && x == 'C' ? 'O' : x
    end

    function ligate(x::Char)
        rand() < s.l_rate && x == 'O' ? 'L' : x
    end 
    

    function cleave(pos::Int,state::Char,sequence::String,aas::String)
        p = rand()
        ### verify w/ oz that nothing should happen if either condition is not true
        if p < s.clv_rate && pos < length(sequence)

            if state == 'C' || state == 'O'
                pos += 1
                state = 'E'
            elseif state == 'L'
                aas *= s.basecall(sequence[pos])
                pos += 1
                state = 'E'

            end
        end 

        (pos,state,sequence,aas)

    end
    state = ligate(oligo(click(x.state)))
    pos,state,sequence,aas = cleave(x.pos,state,x.sequence,x.aas)
    

    
    Peptide(sequence,state,x.oligo,x.ligate,aas,pos)
    
end