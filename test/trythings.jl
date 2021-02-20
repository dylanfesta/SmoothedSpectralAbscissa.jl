using Pkg
Pkg.activate(joinpath(@__DIR__,".."))
using LinearAlgebra
using BenchmarkTools, Test
using Cthulhu
using Plots ; theme(:dark)
using SmoothedSpectralAbscissa ; const SSA = SmoothedSpectralAbscissa

##

function easydyn(x0::Vector{R},mat::Matrix{R},tmax::Real,dt::Real=0.01) where R
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

##

n = 350
idm =diagm(0=>fill(1.0,n))
matrand=randn(n,n) ./ sqrt(n) + 10.0I
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
