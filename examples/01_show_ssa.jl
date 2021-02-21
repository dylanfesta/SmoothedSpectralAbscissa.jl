# # Comparison between SSA and SA

# ### Initialization
using Plots,NamedColors
using LinearAlgebra
using Random
using SmoothedSpectralAbscissa ; const SSA=SmoothedSpectralAbscissa
Random.seed!(0);

# the function below generates a random non-normal matrix
function rand_nonnormal(n::Integer,(ud::Real)=1.01)
  mat = randn(n,n) ./ sqrt(n)
  sh = schur(mat)
  upd = diagm(0=>fill(1.0,n),1=>fill(ud,n-1))
  return sh.vectors*upd*sh.Schur*inv(upd)*sh.vectors'
end;

# ### Example matrix generated  parametrically
n = 100
mat1 = rand_nonnormal(n,0.8)
mat2 = randn(n,n) ./ sqrt(n)
mat(θ) = @. θ*mat1 + (1-θ)*mat2

thetas = range(0.0,1.0;length=100);

# Now we look at the spectral abscissa as a function of the parameter
## #src

# ### Set ϵ for the SSA
ssa_eps_vals = [0.005,0.001,0.0005];

# ### Compute and plot the SA
sas = map(θ->SSA.spectral_abscissa(mat(θ)),thetas)
plt=plot(thetas,sas ; color=:black, linewidth=3,lab="SA",
  xlabel="θ",ylabel="SA value",leg=:top)

# ### Compute and plot the SSA
mycols=cgrad([colorant"DarkGreen",colorant"orange"])
for ϵ in ssa_eps_vals
  ssas = map(θ->SSA.ssa_simple(mat(θ),ϵ),thetas)
  plot!(plt,thetas,ssas ; linewidth=2.5 ,
    lab="SSA $ϵ",palette=mycols)
end
plot!(plt,xlabel="θ",ylabel="SA/SSA value",leg=:top)

#=
This plot shows that the SSA is a smoothed version of the SA, and
coverges to it as as ϵ decreases.
Moreover if we reduce the SSA parametrically, we find (approximately) a good minimum for
the SA, as well.
=#


# Literate.markdown("examples/01_show_ssa.jl","docs/src";documenter=true,repo_root_url="https://dylanfesta.github.io/SmoothedSpectralAbscissa.jl") #src
