using Pkg
Pkg.activate(joinpath(@__DIR__,".."))
using LinearAlgebra
using BenchmarkTools, Test
using Cthulhu
using Plots ; theme(:dark)
using SmoothedSpectralAbscissa ; const SSA = SmoothedSpectralAbscissa

# using OrdinaryDiffEq
using StochasticDiffEq

##

function easydyn(mat::Matrix{R},x0::Vector{R},tmax::Real,dt::Real=0.01) where R
  ts = range(0,tmax;step=dt)
  ret = Matrix{R}(undef,length(x0),length(ts))
  for (k,t) in enumerate(ts)
    ret[:,k]=exp(mat.*t)*x0
  end
  retnrm = mapslices(norm,ret;dims=1)[:]
  return ts,ret,retnrm
end

function lessnormal(mat::Matrix{R},d::R) where R
  n=size(mat,1)
  sh = schur(mat)
  upd = diagm(0=>fill(1.0,n),1=>fill(d,n-1))
  return sh.vectors*upd*sh.Schur*inv(upd)*sh.vectors'
end

function dyn_noise(mat::Matrix{R},
    x0::Vector{R},noiselevel::Real,t_end::Real;
    stepsize::Real=0.05,verbose=false) where R<:Real
  f(du,u,p,t) = let _du = mat*u ; copy!(du,_du) ; end
  σ_f(du,u,p,t) = fill!(du,noiselevel)
  prob = SDEProblem(f,σ_f,x0,(0.,t_end))
  solv =  solve(prob,EM();verbose=verbose,dt=stepsize)
  ret_u = hcat(solv.u...)
  retnrm = mapslices(norm,ret_u;dims=1)[:]
  return solv.t,ret_u,retnrm
end
##

n = 88
idm =diagm(0=>fill(1.0,n))
matrand=randn(n,n) ./ sqrt(n) - 1.1I
matrandln = lessnormal(matrand,1.00001)
alloc=SSA.SSAAlloc(n)

##
SSA.ssa_simple!(matrand,nothing,alloc)
SSA.ssa_simple!(matrandln,nothing,alloc)

##

@btime  SSA.ssa_simple!($matrand,nothing,$alloc)
@btime  SSA.ssa_simple!($matrandln,nothing,$alloc)
println("\n\n")
@btime  SSA.ssa_simple_newton!($matrand,nothing,$alloc)
@btime  SSA.ssa_simple_newton!($matrandln,nothing,$alloc)

##
xstart = randn(n)
t_max=50.0
_step=0.1
mymat=matrandln
times,dynntw= dyn_noise(mymat,xstart,0.1,t_max;stepsize=_step)
times,dynntw_nonoise,_= easydyn(mymat,xstart,t_max,_step)

plot(dynntw[1,:],dynntw[2,:] ; leg=false)
plot!(dynntw_nonoise[1,:],dynntw_nonoise[2,:] ; leg=false)
