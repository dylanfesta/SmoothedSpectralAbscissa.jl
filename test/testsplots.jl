using Pkg ; Pkg.activate(joinpath(@__DIR__,".."))
using SmoothedSpectralAbscissa ; const SSA = SmoothedSpectralAbscissa
using Test
using Plots
using BenchmarkTools
using LinearAlgebra
using Random
using Roots

##

"""
    lyap_simple(A::AbstractMatrix)
returns the naive implementation of the Lyapunov objective function and its
derivative. Can be used for testing, plotting, etc
"""
function lyap_simple(A::AbstractMatrix)
    function  myf(s::Float64)::Float64
        P = lyap(A - UniformScaling(s),one(A))
        tr(P)
    end
    function d_myf(s::Float64)::Float64
        sc,Id = UniformScaling(s), one(A)
        P =  lyap(A - sc,Id)
        Q = lyap(A' - sc,Id)
        -2.0tr(P*Q)
    end
    myf,d_myf
end


"""
    ssa_test(A::Matrix{Float64})

This is the easiest implementation, to make sure it is the same as in the
 ocaml code, and try to benchmark it already.
"""
function ssa_test(A::Matrix{Float64})
    myf,d_myf = lyap_simple(A)
    ssa_eps = SSA.default_eps_ssa(A)
    myg(s)= 1.0 / myf(s) - ssa_eps
    d_myg(s) = - d_myf(s) / myf(s)^2
    sa = SSA.spectral_abscissa(A)
    ssa_start = sa + eps(3sa)
    find_zero( (myg,d_myg),ssa_start, Roots.Newton())
    # find_zero( myg ,ssa_start, Order1())
end


##
function ssa_simple(A::AbstractMatrix, with_gradient::Bool,PQ,
        ssa_eps::Union{Nothing,Float64}=nothing)
    _ssa_eps = something(ssa_eps, SSA.default_eps_ssa(A))
    SSA.PQ_init!(A,PQ)
    _sa = SSA.spectral_abscissa(PQ)
    _start = _sa + 1.05_ssa_eps
    g(s)=SSA._g(s,PQ,_ssa_eps,_sa)
    x = range(_start, 10 ; length=200)
    y = g.(x)
    (x,y)
end
##

# check SA
n = 100
mat = randn(n,n) + UniformScaling(7.4517)
Malloc = SSA.SSAAlloc(n)
@test begin
    sa1 = maximum(real.(eigvals(mat)))
    SSA.PQ_init!(mat,Malloc)
    sa2 = SSA.spectral_abscissa(Malloc)
    sa3 = SSA.spectral_abscissa(mat)
    isapprox(sa1,sa2 ; rtol=1E-4) && isapprox(sa1,sa3 ; rtol=1E-4)
end
# now test SSA against simpler version
ssa1 = ssa_test(mat)

SSA.default_eps_ssa(mat) |> typeof
ssa2 = SSA.ssa_simple(mat)
@test isapprox(ssa1,ssa2 ; rtol=1E-4)
# same , but with a different epsilon
mat = randn(n,n) + UniformScaling(1.456)
eps = 0.31313131
ssa1 = ssa_test(mat; ssa_eps=eps)
ssa2 = SSA.ssa_simple(mat,eps)
@test isapprox(ssa1,ssa2 ; rtol=1E-4)


##
const n = 100
Alloc = SSA.SSAAlloc(n)
A = randn(n,n)

SSA.PQ_init!(A,Alloc)
SSA.spectral_abscissa(A)
SSA.spectral_abscissa(Alloc)

ssa_test(A)
xt,yt = ssa_simple(A,false,Alloc)

plot(xt,yt ; leg = false, linewidth=3)

using Cthulhu

@code_warntype SSA.ssa_simple(A)

descend(SSA.get_P!,
        Tuple{Float64,SSA.SSAAlloc{Array{Float64,2},Array{Complex{Float32},1}}})
