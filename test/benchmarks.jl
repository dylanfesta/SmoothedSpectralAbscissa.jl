
using Pkg
Pkg.activate(joinpath(@__DIR__,".."))
using LinearAlgebra
using BenchmarkTools, Test
using Cthulhu
using Plots ; theme(:dark)
using SmoothedSpectralAbscissa ; const SSA = SmoothedSpectralAbscissa


function rand_gauss(n::Integer)
  return randn(n,n) ./ sqrt(n)
end
function rand_nonnormal(n::Integer,(ud::Real)=1.01)
  mat = rand_gauss(n)
  @show SSA.spectral_abscissa(mat)
  mat = mat -  (1.1*SSA.spectral_abscissa(mat))*I
  sh = schur(mat)
  upd = diagm(0=>fill(1.0,n),1=>fill(ud,n-1))
  return sh.vectors*upd*sh.Schur*inv(upd)*sh.vectors'
end


##

sizes = [10,30,50,100,150,200]

ssa_times = similar(sizes,Float64)
ssa_allocs = similar(ssa_times)


for (k,n) in enumerate(sizes)
  alloc = SSA.SSAAlloc(n)
  A=rand_gauss(n)
  be = @benchmark SSA.ssa!($A,nothing,$alloc)
  ssa_times[k] = median(be.times)
  ssa_allocs[k] = be.allocs
end

##
