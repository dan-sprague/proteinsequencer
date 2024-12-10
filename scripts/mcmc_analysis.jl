using DataFrames,CSV

using Turing
using FillArrays
using Lux
using Plots
using Tracker
using Functors

using LinearAlgebra
using Random


df = CSV.read("../mcmc_sim.csv",DataFrame)


x = Matrix(df[:,[3,6:9...]])
y = df.maxContigLength

rng = Random.default_rng()


@model function poisson_regression(x, y)
    θ ~ MvNormal(zeros(5),1)
    b ~ Normal(0,1)
    m = sum((θ .* x) .+ b )
    λ = exp.(m)
    @. y ~ Poisson(λ)
end;

x = (x .- mean(x;dims=1)) ./ std(x;dims=1)
m = poisson_regression(permutedims(x),y)
chain = sample(m, NUTS(), 2_000) 

df[:,[3,6:9...]]

θ = MCMCChains.group(chain, :θ).value
using StatsPlots
using LaTeXStrings
boxplot(exp.(θ[:,:,1]),side=:left,bandwidth=.15,color=:white,linecolor=:black,size=(150,300),label=:none,ms=1)
ylabel!(L"\exp{θ}")
xticks!(1:5,["BaseCaller","Click","Oligo","Ligate","Cleave"],xrotation=90)
savefig("../poisson_reg_analysis.pdf")
xtest = df[:,[3,6:9...]]

xtest.click .= 0.9
xtest.oligo .= 0.94
xtest.ligate .= 0.85
xtest.cleave .= 0.95

xtest = permutedims(Matrix(xtest))

yy = [map(t -> nn_forward(t,θ[i,:])[1],eachcol(permutedims(x))) for i in 200:2000]
yy = reduce(hcat,yy)


using Zygote


Zygote.gradient((x,θ) -> (y[100] - nn_forward(x,θ)[1])^2,x[100,:],θ[400,:])