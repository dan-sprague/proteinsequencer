using DataFrames,CSV
using Plots,StatsPlots 
using StatsBase

df = CSV.read("./digest.csv",DataFrame)
df.proportionOfReads .= round.(df.proportionOfReads;digits=4)
df.maxContigPercIden ./= 100

meanAccCov = combine(groupby(df,[:acc,:proportionOfReads]), :fractionMaxContigCoverage => mean,:fractionMaxContigCoverage => std)
sort!(meanAccCov,:acc)

gdf = groupby(meanAccCov,:proportionOfReads)


pl = []

for g in gdf
    sort!(g,:acc)
    p = @df g plot(:acc,:fractionMaxContigCoverage_mean,yerr=:fractionMaxContigCoverage_std,
    title=first(g.proportionOfReads),label=:none,titlefont=(10))
    xlims!(0,1)
    push!(pl,p)
end

plot(pl...)
savefig("./acc_propofref.pdf")


meanAccCov = combine(groupby(df,[:acc,:proportionOfReads]), :maxContigPercIden => mean,:maxContigPercIden => std)
sort!(meanAccCov,:acc)

gdf = groupby(meanAccCov,:proportionOfReads)


pl = []

for g in gdf
    sort!(g,:acc)
    p = @df g plot(:acc,:maxContigPercIden_mean,yerr=:maxContigPercIden_std,
    title=first(g.proportionOfReads),label=:none,titlefont=(10))
    xlims!(0,1)
    ylims!(0,1)
    push!(pl,p)
end

plot(pl...)
savefig("./acc_percid.pdf")


### perfect sequencing

df = CSV.read("./perfect_sequencing_normalcycle.csv",DataFrame)
df.proportionOfReads .= round.(df.proportionOfReads;digits=4)
df.maxContigPercIden ./= 100

using GaussianProcesses
scatter(df.cycles,df.fractionMaxContigCoverage)

using StatsBase


meanAccCov = combine(groupby(df,[:cycles,:proportionOfReads]), :fractionMaxContigCoverage => mean,:fractionMaxContigCoverage => std)

x = permutedims(Matrix(meanAccCov[:,[1,2]]))
y = meanAccCov[:,3]

gdf = groupby(meanAccCov,:proportionOfReads)

pl = []

for g in gdf
    sort!(g,:cycles)
    p = @df g plot(:cycles,:fractionMaxContigCoverage_mean,yerr=:fractionMaxContigCoverage_std,
    title=first(g.proportionOfReads),label=:none,titlefont=(12))
    xticks!([6,10,15,30])
    xlims!(0,35)
    push!(pl,p)
end

plot(pl...)

savefig("./perfectsequencingnormal.pdf")


meanAccId = combine(groupby(df,[:cycles,:proportionOfReads]), :maxContigPercIden => mean,:maxContigPercIden => std)

x = permutedims(Matrix(meanAccId[:,[1,2]]))
y = meanAccId[:,3]

gdf = groupby(meanAccId,:proportionOfReads)

pl = []

for g in gdf
    sort!(g,:cycles)
    p = @df g plot(:cycles,:maxContigPercIden_mean,yerr=:maxContigPercIden_std ./ sqrt(5),
    title=first(g.proportionOfReads),label=:none,titlefont=(12))
    xticks!([6,10,15,30])
    xlims!(0,35)
    ylims!(0,1)
    push!(pl,p)
end

plot(pl...)

savefig("./perfectsequencingnormal_percid.pdf")




df = CSV.read("./perfect_sequencing_longcycle.csv",DataFrame)
df.proportionOfReads .= round.(df.proportionOfReads;digits=4)

using GaussianProcesses
scatter(df.cycles,df.fractionMaxContigCoverage)

using StatsBase


meanAccCov = combine(groupby(df,[:cycles,:proportionOfReads]), :fractionMaxContigCoverage => mean,:fractionMaxContigCoverage => std)

x = permutedims(Matrix(meanAccCov[:,[1,2]]))
y = meanAccCov[:,3]

gdf = groupby(meanAccCov,:proportionOfReads)

pl = []

for g in gdf
    sort!(g,:cycles)
    p = @df g plot(:cycles,:fractionMaxContigCoverage_mean,yerr=:fractionMaxContigCoverage_std,
    title=first(g.proportionOfReads),label=:none,titlefont=(12))
    xticks!([30,50,100])
    xlims!(0,105)
    ylims!(0,2)
    push!(pl,p)
end

plot(pl...)

savefig("./perfectsequencinglong.pdf")


meanAccId = combine(groupby(df,[:cycles,:proportionOfReads]), :maxContigPercIden => mean,:maxContigPercIden => std)

x = permutedims(Matrix(meanAccId[:,[1,2]]))
y = meanAccId[:,3]

gdf = groupby(meanAccId,:proportionOfReads)

pl = []

for g in gdf
    sort!(g,:cycles)
    p = @df g plot(:cycles,:maxContigPercIden_mean,yerr=:maxContigPercIden_std ./ sqrt(5),
    title=first(g.proportionOfReads),label=:none,titlefont=(12))
    xticks!([30,50,100])
    xlims!(0,105)
    ylims!(0,1)
    push!(pl,p)
end

plot(pl...)

savefig("./perfectsequencinglong_percid.pdf")


### CURRENT PARAMS ####
using Plots
using Distributions
using StatsPlots 
using DataFrames,CSV

types = Dict(
    :ratio => String,
)
df = DataFrame(CSV.File("./digest.csv";types))
df.proportionOfReads .= round.(df.proportionOfReads;digits=4)



gdf = combine(groupby(df,[:cycles,:ratio,:nProts,:proportionOfReads,:acc]), :maxContigLength => mean,:maxContigLength => std)
sort!(gdf,:acc)
gdf = groupby(gdf,:acc)

pl = []
prop = []
#g = plot()
for grp in gdf
    sort!(grp,:proportionOfReads)
    p = scatter(1:size(grp,1),grp.maxContigLength_mean,yerr=grp.maxContigLength_std / sqrt(5),
            group=grp.cycles,title="acc = $(first(grp.acc))",bottommargin=2Plots.mm
    )

    r = string.(length.(split.(grp.ratio,'.')))

    ylims!(0,1000)
    xticks!(1:size(grp,1), "n =" .* r .* " | " .* string.(grp.proportionOfReads),xrotation=45)
    push!(pl,p)
end

plot(pl...,size=(200,300),ytickfont=(10),xtickfont=(10),xrotation=90,bottommargin=10Plots.mm)
savefig("./maxcontiglength_ecoli.pdf")

gdf = combine(groupby(df,[:cycles,:ratio,:proportionOfReads,:acc]), :maxContigLength => mean,:maxContigLength => std,:maxContigPercIden => mean,:maxContigPercIden => std,
:maxBlastTaxIdFracCorrect => mean, :maxBlastTaxIdFracCorrect => std)
gdf = gdf[gdf.acc .>= 0.8 .&& gdf.cycles .>= 15,:]
sort!(gdf,[:cycles,:acc])
gdf = groupby(gdf,[:cycles,:acc])
pl = []
g = plot()
for k in keys(gdf)
    grp = gdf[k]
    idx = grp.maxContigLength_mean .== maximum(grp.maxContigLength_mean)
    μ= grp.maxContigLength_mean[idx][1]
    σ = grp.maxContigLength_std[idx][1] ./ sqrt(5)
    xx = 1:2000
    xxx = []
    for i in 1:2000
        r = rand(Normal(μ[1],max(20,σ)))
        r = rand(Poisson(max(1,r)))
        push!(xxx,[sum(1:x .< r ) / x for x in 1:2000])
    end
    xxx = reduce(hcat,xxx)
    errorline!(xxx,xlabel="Reference Length (aa)",ylabel="Expected Contiguous\nCoverage",
        label="cycles = $(k.cycles)\nbasecaller = $(k.acc)\n% Identity = $(round(grp[idx,:maxContigPercIden_mean][1];digits=2)) ± $(round(grp[idx,:maxContigPercIden_std][1] / sqrt(5);digits=2))\n% Correct TaxID = $(round(grp[idx,:maxBlastTaxIdFracCorrect_mean][1];digits=2) * 100) ± $(round(grp[idx,:maxBlastTaxIdFracCorrect_std][1] / sqrt(5);digits=2))",
        linewidth=2.0)
    ylims!(0,1)
    push!(pl,g)
end


plot(g,size=(500,300),leftmargin=-15Plots.mm,tickfont=(8),legendfont=(6),
guidefont=(9),titlefont=(10),xrotation=45,legend=:outerleft,bottommargin=2Plots.mm)


savefig("acc_proportion.pdf")

sort!(gdf,[:acc,:nProts])
g = groupby(gdf,[:acc,:nProts])

pl = []
for grp in g
    acc = first(grp.acc)
    nprots = first(grp.nProts)
    ylabel = nprots == 2 ? "acc = $acc\nMax Length" : ""
    xlabel = acc == 0.95 ? "Read Proportion" : ""
    sort!(grp,:proportionOfReads)
    push!(pl,plot(grp.proportionOfReads,grp.maxContigLength_mean,group=grp.cycles,ribbon=grp.maxContigLength_std ./ sqrt(5),
    linewidth=1.5,ylabel=ylabel,xlabel=xlabel),)
end

title1 =  plot(title = "Proteins = 2", grid = false, showaxis = false, bottom_margin = -150Plots.px)
title2 =  plot(title = "Proteins = 3", grid = false, showaxis = false, bottom_margin = -150Plots.px)
plot(title1,title2,pl...,layout=grid(6,2),size=(350,800),xrotation=45,leftmargin=5Plots.mm,ylims=(0,500),
legend=:topleft,tickfont=(8),legendfont=(6),
guidefont=(9),titlefont=(10))
savefig("FIGURE_5.pdf")

gdf.maxContigLength_sem = gdf.maxContigLength_std ./ sqrt(5)


CSV.write("TABLE_5.csv",gdf)