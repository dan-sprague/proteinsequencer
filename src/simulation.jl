"""
    simulate!(sequencer::Sequencer,peptides,X;NCYCLES = 15)

Pass the `peptides` through the in silico `sequencer` for `NCYCLES`. Saves summary statistics to arrays provided in `X`. Modifies `peptides` and `X` in place. 
"""
function simulate!(sequencer::Sequencer,peptides,X;NCYCLES=15)

    μ_length,σ_length,μ_pos,σ_pos,μ_deletions,σ_deletions,dels = X

    @inbounds for i ∈ 1:NCYCLES
        @. peptides = sequencer(peptides)
        l = [length(x.aas) for x in peptides]
        pos = [x.pos for x in peptides]
        @. dels = pos - l

        μ_length[i] = mean(l)
		σ_length[i] = std(l)

		μ_pos[i] = mean(pos)
		σ_pos[i] = std(pos)


		μ_deletions[i] = mean(dels)
		σ_deletions[i] = std(dels)
    end
end 

"""
    plot(X)

Plots summary statistics for read length, position, and deletions.
"""
function plot(X)
    μ_length,σ_length,μ_pos,σ_pos,μ_deletions,σ_deletions,dels = X
    
    p1 = Plots.plot(μ_length,ribbon=σ_length,label=:none)
	ylabel!("AAs sequenced")

	p2 = Plots.plot(μ_pos,ribbon=σ_pos,label=:none)
	ylabel!("Position in Sequence")

	p3 = Plots.plot(μ_deletions,ribbon=σ_deletions,label=:none)
	xlabel!("Cycles")
	ylabel!("Deletions")


	layout = @layout [a;b;d]

	g = Plots.plot(p1,p2,p3,
        layout=layout,
        linewidth=2,
        fillalpha=0.3,
        tickfont=(10,"Computer Modern"),
        guidefont=(12,"Computer Modern"),
        size=(300,700),
        leftmargin=10Plots.mm
    )

    savefig(g,"./test.pdf")
end
