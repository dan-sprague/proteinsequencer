function pass_report(path)
    files = glob("$path/reads.fasta.*")
    p = r"pid\d+\."
    file_groups = Dict{String, Vector{String}}()
	for file in files
        m = match(p, file)
        if m !== nothing
            key = m.match
            if haskey(file_groups, key)
                push!(file_groups[key], file)
            else
                file_groups[key] = [file]
            end
        else
            continue
        end
    
        file_groups[key] = unique(file_groups[key])
    end

	statisticsDict = Dict(zip(collect(keys(file_groups)),copy.(repeat([Dict(
		"N" => 0.0,
		"short" => 0.0,
		"singlets" => 0.0,
		"efficiency" => 0.0,
		"contigs" => 0.0,
		"maxlength" => 0,
		"meanlength" => 0.0,
		"meancoverage" => 0.0,
		)],length(keys(file_groups))))
	)
	)
	pattern = "total sequences"
	for assembly in keys(file_groups)
		
		

		open(file_groups[assembly][2]) do file
			for line in eachline(file)
				if occursin(pattern,line)
					statisticsDict[assembly]["N"] = parse(Int,split(line)[1])
				end
			end
		end

		open(file_groups[assembly][3]) do file

			statisticsDict[assembly]["short"] =  countlines(file) / statisticsDict[assembly]["N"]
		end

		open(file_groups[assembly][4]) do file
			statisticsDict[assembly]["singlets"] = (countlines(file) ÷ 2) / statisticsDict[assembly]["N"]

		end

		open(file_groups[assembly][1]) do file
			C = 0
			maxL = 1
			L = 0.0
			for line in eachline(file)
				if startswith(line,">")
					statisticsDict[assembly]["meancoverage"] += parse(Float32,split(line,"|cov")[end])
                    C += 1 # fasta file... every other line is a count :)
				else
					l = length(line)
					l > maxL ? maxL = l : maxL = maxL
					L += l
				end
			end

			statisticsDict[assembly]["contigs"] = C
			statisticsDict[assembly]["meancoverage"] /= C
			statisticsDict[assembly]["meanlength"]	= round(L / C,digits = 3)
			statisticsDict[assembly]["maxlength"] = maxL
			statisticsDict[assembly]["efficiency"] = (statisticsDict[assembly]["meancoverage"] * statisticsDict[assembly]["contigs"]) / statisticsDict[assembly]["N"]

            statisticsDict[assembly]["efficiency"] = round(statisticsDict[assembly]["efficiency"],digits = 3)
            statisticsDict[assembly]["meancoverage"] = round(statisticsDict[assembly]["meancoverage"],digits = 3)
            statisticsDict[assembly]["short"] = round(statisticsDict[assembly]["short"],digits=3)
            statisticsDict[assembly]["singlets"] = round(statisticsDict[assembly]["singlets"],digits=3)
		end

	end

	flattened_data = []
    for (group, files) in statisticsDict
        for (file, value) in files
            push!(flattened_data,[group, file,value])
        end
    end

    df = DataFrame(permutedims(reduce(hcat,flattened_data),(2,1)),:auto)

    μ_maxl = mean(df.x3[df.x2 .== "maxlength"])
    stderr_maxl = std(df.x3[df.x2 .== "maxlength"]) / sqrt(length(df.x3[df.x2 .== "maxlength"]))

    @show "Mean Max Contig Length ± Standard Error: $(round(μ_maxl,digits=2)) ± $(round(stderr_maxl,digits=4))"


    μ_eff = mean(df.x3[df.x2 .== "efficiency"])
    stderr_eff = std(df.x3[df.x2 .== "efficiency"]) / sqrt(length(df.x3[df.x2 .== "efficiency"]))

    @show "Mean Efficiency ± Standard Error: $(round(μ_eff,digits=2)) ± $(round(stderr_eff,digits=4))"

	df 
end