using SmoothedSpectralAbscissa ; const SSA=SmoothedSpectralAbscissa
using LinearAlgebra, Random
using Roots,Calculus
using Test

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

This is the easiest implementation, to make sure it is the same as in the
 ocaml code, and try to benchmark it already.
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
#
# function _dg(s,PQ::SSA.SSAAlloc)
#     SSA._get_P!(s,PQ)
#     SSA._get_Q!(s,PQ)
#     fs = tr(PQ.P)
#     mdfs = 2.0*_trace_prod(PQ.P,PQ.Q)
#     return mdfs / fs^2
# end
# function get_ssa_froot(A::AbstractMatrix, ssa_eps::Float64=0.0)
#     _ssa_eps = ssa_eps > 0.0 ? ssa_eps : 0.01 * 150.0 / size(A,1) #scales as 1/n
#     PQ=SSA.SSAAlloc(size(A,1))
#     SSA.PQ_prepare(A,PQ)
#     g(s)=_g(s,PQ, _ssa_eps)
#     dg(s) = _dg(s,PQ)
#     return (g,dg)
# end
#
# relerror(a,b) = let apb=a+b ; apb == 0.0 ? 0.0 : 2*abs(a-b)/(apb) ; end
#
# function gradient_test(A::AbstractMatrix,s::AbstractVector, ssa_eps::Float64=0.0)
#     g,dg = get_ssa_froot(A,ssa_eps)
#     an = dg.(s)
#     num = map( ss -> Calculus.gradient(g,ss) , s )
#     err = broadcast(relerror, an,num)
#     an,num,err
# end

# function f_for_gradient(vmat)
#     n = sqrt(length(vmat))
#     @assert isinteger(n)
#     n = Int64(n)
#     mat = copy(reshape(vmat,(n,n)))
#     return ssa_test(mat)
# end
function f_for_gradient(vmat)
    n = sqrt(length(vmat))
    @assert isinteger(n)
    n = Int64(n)
    mat = copy(reshape(vmat,(n,n)))
    return SSA.ssa_simple(mat)
end

#

@testset "SA and SSA" begin
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
    @test begin
        ssa1 = ssa_test(mat)
        ssa2 = SSA.ssa_simple(mat)
        isapprox(ssa1,ssa2 ; rtol=1E-4)
    end
    # same , but with a different epsilon
    mat = randn(n,n) + UniformScaling(1.456)
    eps = 0.31313131
    @test begin
        ssa1 = ssa_test(mat; ssa_eps=eps)
        ssa2 = SSA.ssa_simple(mat,eps)
        isapprox(ssa1,ssa2 ; rtol=1E-4)
    end
end

@testset "SSA gradient" begin
    n = 26
    mat = rand(n,n) + UniformScaling(2.222)
    @test begin
        ssa1 = ssa_test(mat)
        ssa2 = SSA.ssa_simple(mat)
        isapprox(ssa1,ssa2 ; rtol=1E-4)
    end
    @test begin
        ssa,grad_an = SSA.ssa_simple_withgradient(mat)
        grad_num = Calculus.gradient(f_for_gradient,mat[:])
        all(isapprox.(grad_an[:],grad_num ; rtol=1E-3) )
    end
end
