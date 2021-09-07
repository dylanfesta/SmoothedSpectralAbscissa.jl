
using Test
using Plots
using BenchmarkTools
using LinearAlgebra
using Random
using Roots

push!(LOAD_PATH, abspath(@__DIR__,".."))

using SmoothedSpectralAbscissa ; const SSA = SmoothedSpectralAbscissa

##


"""
    lyap_simple!(A::AbstractMatrix)
returns the naive implementation of the Lyapunov objective function and its
derivative. Can be used for testing, plotting, etc
"""
function lyap_simple!(A::AbstractMatrix)
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

This is the easiest implementation
for testing and benchmarking.
"""
function ssa_test(A::Matrix{Float64} ; ssa_eps=nothing)
    myf,d_myf = lyap_simple!(A)
    ssa_eps = something(ssa_eps,SSA.default_eps_ssa(A))
    myg(s)= 1.0 / myf(s) - ssa_eps
    d_myg(s) = - d_myf(s) / myf(s)^2
    sa = SSA.spectral_abscissa(A)
    ssa_start = sa + eps(3sa)
    find_zero( (myg,d_myg),ssa_start, Roots.Newton())
end

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
