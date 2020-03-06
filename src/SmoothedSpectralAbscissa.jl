module SmoothedSpectralAbscissa
using LinearAlgebra, Roots


spectral_abscissa(A::AbstractMatrix) =  maximum(real.(eigvals(A)))

# trace of matrix product
function trace_of_product(A::AbstractMatrix,B::AbstractMatrix)
    n=size(A,1)
    acc=0.0
    for i in 1:n, j in 1:n
        @inbounds @fastmath acc += A[i,j]*B[j,i]
    end
    return acc
end

function subtract_diagonal!(A,s)
    n=size(A,1)
    @simd for i in 1:n
        @inbounds A[i,i] -= s
    end
    return nothing
end

"""
        default_eps_ssa(A::AbstractMatrix{Float64}) -> eps_ssa

default value for the ``\\varepsilon`` used to compute the SSA
"""
default_eps_ssa(A::AbstractMatrix) = 0.01 * 150.0 / size(A,1)


struct SSAAlloc{M<:AbstractMatrix,V<:AbstractVector}
    R::M
    R_alloc::M
    Rt::M
    Z::M
    Zt::M
    D::M
    D_alloc::M
    Dt::M
    At::M
    P::M
    Q::M
    YZt_alloc::M
    G::M
    A_eigvals::V
end
"""
        SSAAlloc{M<:AbstractMatrix,V<:AbstractVector}(n::Integer) -> SSAAlloc

Allocation structure that makes SSA computation mode efficient. `n` is the size of
the matrix (or matrices) on which one wants to compute the SSA
"""
function SSAAlloc(n::Integer)
    return SSAAlloc(
       map( _ -> Matrix{Float64}(undef,n,n),1:13)... , Vector{ComplexF32}(undef,n))
end

"""
        PQ_init!(A::AbstractMatrix,PQ::SSAAlloc) -> nothing
Prepares SSAAlloc with the matrix on which one wants to compute the SSA.
Schur decompositions are performed here.
"""
function PQ_init!(A::AbstractMatrix,PQ::SSAAlloc)
    At=transpose!(PQ.At,A)
    F = schur(A)
    copyto!(PQ.R,F.T)
    copyto!(PQ.Z,F.Z)
    copyto!(PQ.A_eigvals,F.values)
    F = schur!(At)
    copyto!(PQ.Rt,F.T)
    copyto!(PQ.Zt,F.Z)
    mul!(PQ.D,transpose(PQ.Z),PQ.Z,-1.0,0.0)
    mul!(PQ.Dt,transpose(PQ.Zt),PQ.Zt,-1.0,0.0)
    return nothing
end

function spectral_abscissa(PQ::SSAAlloc)
    maximum(real.(PQ.A_eigvals))
end
function get_P_or_Q!(PorQ,s,Z,R,D,YZt_alloc)
    subtract_diagonal!(R,s)
    Y, scale = LAPACK.trsyl!('N','T', R, R, D)
    mul!(YZt_alloc,Y,transpose(Z))
    mul!(PorQ,Z,YZt_alloc,inv(scale),0.0)
    return nothing
end
function get_P!(s, PQ::SSAAlloc)
    R=copyto!(PQ.R_alloc,PQ.R)
    D=copyto!(PQ.D_alloc,PQ.D)
    get_P_or_Q!(PQ.P,s,PQ.Z,R,D,PQ.YZt_alloc)
    return nothing
end
function get_Q!(s, PQ::SSAAlloc)
    R=copyto!(PQ.R_alloc,PQ.Rt)
    D=copyto!(PQ.D_alloc,PQ.Dt)
    get_P_or_Q!(PQ.Q,s,PQ.Zt,R,D,PQ.YZt_alloc)
    return nothing
end

function ssa_simple_obj(s,PQ::SSAAlloc ,ssa_eps,sa)
    get_P!( max(sa,s) ,PQ) # SSA is not defined below SA !
    return inv(tr(PQ.P)) - ssa_eps
end


"""
        ssa_simple!(A::AbstractMatrix, with_gradient::Bool,PQ::SSAAlloc ,
            ssa_eps::Union{Nothing,Float64}=nothing) -> ssa::Float64

Computes the Smoothed Spectral Abscissa (SSA) of matrix A in its simplest formulation.
For efficiency, `A` is reallocated in the computation. If `ssa_eps` is not specified, a default value is used.
"""
function ssa_simple!(A::AbstractMatrix, with_gradient::Bool,PQ::SSAAlloc ,
        ssa_eps::Union{Nothing,Float64}=nothing)
    _ssa_eps = something(ssa_eps, default_eps_ssa(A))
    PQ_init!(A,PQ)
    _sa = spectral_abscissa(PQ)
    _start = _sa + 1.05_ssa_eps
    g(s)=ssa_simple_obj(s,PQ,_ssa_eps,_sa)
    s_star::Float64 = find_zero( g, _start , Order2())
    if with_gradient
        get_P!(s_star,PQ)
        get_Q!(s_star,PQ)
        cc=inv(trace_of_product(PQ.Q,PQ.P))
        mul!(PQ.G,PQ.Q,PQ.P,cc,0.0)
    end
    return s_star
end
# short versions for 1-shot computation and testing
"""
        ssa_simple(A::AbstractMatrix , ssa_eps::Union{Nothing,Float64}=nothing)
            -> ssa::Float64

Computes the Smoothed Spectral Abscissa (SSA) of matrix A in its simplest formulation
If `ssa_eps` is not specified, a default value is used.

In case you need to call the SSA iteratively, consider the version that acts on
pre-allocated memory.
"""
function ssa_simple(A::AbstractMatrix , ssa_eps::Union{Nothing,Float64}=nothing)
    PQ=SSAAlloc(size(A,1))
    _epsssa = something(ssa_eps, default_eps_ssa(A))
    return ssa_simple!(copy(A),false,PQ , _epsssa)
end

"""
      ssa_simple_withgradient(A::AbstractMatrix , ssa_eps::Union{Nothing,Float64}=nothing)
            -> (ssa,gradient_matrix)::Tuple{Float64,Matrix}
Computes the Smoothed Spectral Abscissa (SSA) of matrix A in its simplest formulation,
and also its gradient, for each matrix element.
If `ssa_eps` is not specified, a default value is used.

In case you need to call the SSA iteratively (e.g. gradient-based optimization),
 consider the version that acts on pre-allocated memory.
"""
function ssa_simple_withgradient(A::AbstractMatrix, ssa_eps::Union{Nothing,Float64}=nothing)
    PQ=SSAAlloc(size(A,1))
    epsssa = something(ssa_eps, default_eps_ssa(A))
    ssa = ssa_simple!(copy(A),true,PQ,epsssa)
    return ssa,PQ.G
end

mutable struct SSAFun{F<:Real}
    s_old::F
    tr_p::F
end
SSAFun() = SSAFun(-1.0,-1.0)
function __g(s,PQ::SSAAlloc ,stuff::SSAFun, ssa_eps,sa)
    s = max(sa,s) # not defined below SA !
    _get_P!(s,PQ)
    tr_p = tr(PQ.P)
    stuff.tr_p=tr_p
    stuff.s_old=s
    return inv(tr_p) - ssa_eps
end
function __dg(s,PQ::SSAAlloc,stuff::SSAFun,sa)
    s = max(sa,s) # not defined below SA !
    do_stuff =  s != stuff.s_old
    do_stuff && _get_P!(s,PQ)
    fs = do_stuff ?  tr(PQ.P) : stuff.tr_p
    _get_Q!(s,PQ)
    mdfs = 2.0*_trace_prod(PQ.P,PQ.Q)
    return mdfs / fs^2
end
"""
This uses the Newton method and the analytic derivative ,
seems as fast as the Order16 method
"""
function ssa_simple_newton(A::AbstractMatrix,
        with_gradient::Bool,
        PQ::SSAAlloc , ssa_eps::Union{Nothing,Float64}=nothing)
    _ssa_eps = something(ssa_eps, default_eps_ssa(A))
    PQ_init!(A,PQ)
    _sa = spectral_abscissa(PQ)
    _start = _sa + 1.05_ssa_eps
    stuff=SSAFun()
    g(s)=__g(s,PQ,stuff, _ssa_eps,_sa)
    dg(s) = __dg(s,PQ,stuff,_sa)
    s_star = find_zero( (g,dg),_start, Roots.Newton() )
    if with_gradient
        P_star = _get_P!(s_star,PQ)
        Q_star= _get_Q!(s_star,PQ)
        cc=inv(_trace_prod(Q_star,P_star))
        mul!(PQ.G,Q_star,P_star,cc,0.0)
    end
    return s_star
end

function ssa_simple_newton(A::AbstractArray , ssa_eps::Float64=0.0)
    PQ=SSAAlloc(size(A,1))
    return ssa_simple_newton(A,false,PQ , ssa_eps)
end



end # module
