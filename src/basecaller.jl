"""
    Basecaller(acc::Float32,AA::Vector{Char})

Defines a basecalling simulator based on a global accuracy `acc` for the list of amino acids `AA`. This can be extended for more nuanced error rates in the future.

    ```julia
    basecaller = BaseCaller(0.8,AAs)

    for aa in peptide
        basecall(aa)
    end
    ```
"""
struct BaseCaller
    acc::Float32
    AA::Vector{Char}

    function BaseCaller(acc,AA)
        @assert acc <= 1 && acc >= 0
        new(acc,AA)
    end
end

function (b::BaseCaller)(aa)
    rand() < b.acc ? aa : sample(b.AA[b.AA .!= aa])
end

