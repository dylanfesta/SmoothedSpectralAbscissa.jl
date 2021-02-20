module SmoothedSpectralAbscissa
using LinearAlgebra, Roots


@inline function spectral_abscissa(A::Matrix{R}) where R
   convert(R,maximum(real.(eigvals(A))))
end

# trace of matrix product
function trace_of_product(A::Matrix{R},B::Matrix{R}) where R
    n=size(A,1)
    acc=0.0
    for i in 1:n, j in 1:n
        @inbounds @fastmath acc += A[i,j]*B[j,i]
    end
    return acc
end

function subtract_diagonal!(A::Matrix{R},s::R) where R
    n=size(A,1)
    @simd for i in 1:n
        @inbounds A[i,i] -= s
    end
    return nothing
end

"""
        default_eps_ssa(A::Matrix{Float64}) -> eps_ssa

default value for the ``\\varepsilon`` used to compute the SSA. It scales with the matrix
size, so that it takes the value 0.05 for a ``150 \\times 150`` matrix.
"""
default_eps_ssa(A::Matrix{<:Real}) = 0.01 * 150.0 / size(A,1)


struct SSAAlloc{R,C}
    R::Matrix{R}
    R_alloc::Matrix{R}
    Rt::Matrix{R}
    Z::Matrix{R}
    Zt::Matrix{R}
    D::Matrix{R}
    D_alloc::Matrix{R}
    Dt::Matrix{R}
    At::Matrix{R}
    P::Matrix{R}
    Q::Matrix{R}
    YZt_alloc::Matrix{R}
    A_eigvals::Vector{C}
end
"""
        SSAAlloc(n::Integer) -> SSAAlloc

Memory allocation to avoid re-allocating space for repeated SSA computations.
# Arguments
- `n::Integer`: the size of the matrix (or matrices) on which the SSA will be computed
"""
function SSAAlloc(n::Integer)
  return SSAAlloc(
   map( _ -> Matrix{Float64}(undef,n,n),1:12)... , Vector{ComplexF32}(undef,n))
end

"""
        SSAAlloc(A::Matrix) -> SSAAlloc

Memory allocation to avoid re-allocating space for repeated SSA computations.
This corresponds to `SSAAlloc(size(A,1))`
"""
function SSAAlloc(A::Matrix{<:Real})
  return SSAAlloc(size(A,1))
end

"""
        PQ_init!(PQ::SSAAlloc,A::Matrix) -> nothing
Stores some (costly) pre-computed quantities on SSAAlloc.
Should be called only once before the root-finding procedure.
Schur decompositions are performed here.
"""
function PQ_init!(PQ::SSAAlloc{R,C},A::Matrix{R}) where {R,C}
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

@inline function spectral_abscissa(PQ::SSAAlloc{R,C}) where {R,C}
    return convert(R,maximum(real.(PQ.A_eigvals)))
end
function get_P_or_Q!(PorQ::M,s::Real,Z::M,R::M,D::M,
        YZt_alloc::M) where M<:Matrix{<:Real}
    subtract_diagonal!(R,s)
    Y, scale = LAPACK.trsyl!('N','T', R, R, D)
    mul!(YZt_alloc,Y,transpose(Z))
    mul!(PorQ,Z,YZt_alloc,inv(scale),0.0)
    return nothing
end
function get_P!(s::T, PQ::SSAAlloc{T,C}) where {T,C}
    R=copyto!(PQ.R_alloc,PQ.R)
    D=copyto!(PQ.D_alloc,PQ.D)
    get_P_or_Q!(PQ.P,s,PQ.Z,R,D,PQ.YZt_alloc)
    return nothing
end
function get_Q!(s::T, PQ::SSAAlloc{T,C}) where {T,C}
    R=copyto!(PQ.R_alloc,PQ.Rt)
    D=copyto!(PQ.D_alloc,PQ.Dt)
    get_P_or_Q!(PQ.Q,s,PQ.Zt,R,D,PQ.YZt_alloc)
    return nothing
end

function ssa_simple_obj(s::R,PQ::SSAAlloc{R,C},ssa_eps::R,sa::R) where {R,C}
    get_P!( max(sa,s) ,PQ) # SSA is not defined below SA !
    return inv(tr(PQ.P)) - ssa_eps
end


"""
    ssa_simple!(A,grad,PQ,ssa_eps=nothing) -> ssa

Computes the smoothed spectral abscissa (SSA) of matrix A (with identity input-output weighting matrices) and its gradient (optionally).
For efficiency, `A` is reallocated in the computation.
If `ssa_eps` is not specified, the default value is set by `default_eps_ssa(A)`

# Arguments
- `A::Matrix{<:Real}`: A square matrix
- `grad::Union{Nothing,Matrix{<:Real}}`: `nothing`  skips the gradient computation. To compute the gradient, provide a matrix `grad` of the same size as A. The matrix will be rewritten in-place with entries ``\\frac{\\partial\\; \\tilde{\\alpha}_\\varepsilon(A)}{\\partial\\; A_{i,j}}``
- `PQ::SSAAlloc` : Pre-allocated matrices for iterative computations. Can be generated by `SSAAlloc(A)`
- `ssa_eps::Union{Nothing,R}=nothing` The ``\\varepsilon`` parameter associated with the SSA computation. If not specified, it is set by `default_eps_ssa(A)`
# Returns
- `ssa<:Real`: this is ``\\tilde{\\alpha}_\\varepsilon(A)``
"""
function ssa_simple!(A::Matrix{R},grad::Union{Nothing,Matrix{R}},
        PQ::SSAAlloc{R,C},ssa_eps::Union{Nothing,R}=nothing) where {R,C}
    _ssa_eps = something(ssa_eps, default_eps_ssa(A))
    PQ_init!(PQ,A)
    _sa = spectral_abscissa(PQ)
    # _start = _sa + 0.1abs(_sa)
    _start = _sa + 0.5_ssa_eps
    objfun(s)=ssa_simple_obj(s,PQ,_ssa_eps,_sa)
    s_star::Float64 = find_zero(objfun, _start , Order2();maxevals=1_000)
    isnothing(grad) && return s_star
    get_P!(s_star,PQ)
    get_Q!(s_star,PQ)
    cc=inv(trace_of_product(PQ.Q,PQ.P))
    mul!(grad,PQ.Q,PQ.P,cc,0.0)
    return s_star
end


# short versions that also allocates the memory
"""
    ssa_simple(A,ssa_eps=nothing) -> ssa

Computes the smoothed spectral abscissa (SSA) of matrix A (with identity input-output weighting matrices). If `ssa_eps` is not specified, the default value is set by `default_eps_ssa(A)`

# Arguments
- `A::Matrix{<:Real}`: A square matrix
- `ssa_eps::Union{Nothing,R}=nothing` The ``\\varepsilon`` parameter associated with the SSA computation. If not specified, it is set by `default_eps_ssa(A)`
# Returns
- `ssa<:Real`: this is the ``\\tilde{\\alpha}_\\varepsilon(A)``
"""
function ssa_simple(A::Matrix{R},ssa_eps::Union{Nothing,R}=nothing) where R
    PQ=SSAAlloc(size(A,1))
    _epsssa = something(ssa_eps, default_eps_ssa(A))
    return ssa_simple!(copy(A),nothing,PQ , _epsssa)
end

"""
      ssa_simple_withgradient(A,ssa_eps) -> (ssa,gradient_matrix)
Computes the Smoothed Spectral Abscissa (SSA) of matrix A
and also its gradient with respect to each element of A.
If `ssa_eps` is not specified, the default value is set by `default_eps_ssa(A)`

In case you need to call the SSA iteratively (e.g. gradient-based optimization),
please consider pre-allocating memory and using `ssa_simple!(...)`.
"""
function ssa_simple_withgradient(A::Matrix{R},
        ssa_eps::Union{Nothing,R}=nothing) where R
    PQ=SSAAlloc(size(A,1))
    gradmat=similar(A)
    epsssa = something(ssa_eps, default_eps_ssa(A))
    ssa = ssa_simple!(copy(A),gradmat,PQ,epsssa)
    return ssa,gradmat
end


# Find zero with Newton!

function _trace_prod(A::Matrix{R},B::Matrix{R}) where R
  n=size(A,1)
  ret=0.0
  for i in 1:n, j in 1:n
    @inbounds ret += A[i,j]*B[j,i]
  end
  return ret
end


function ssa_simple_obj_newton(s::R,PQ::SSAAlloc{R,C},ssa_eps::R,sa::R) where {R,C}
  _s = max(sa,s)
  get_P!(_s,PQ) # SSA is not defined below SA !
  get_Q!(_s,PQ)
  fs = tr(PQ.P)
  mdfs = 2.0*_trace_prod(PQ.P,PQ.Q)
  obj = inv(fs) - ssa_eps
  dobj = mdfs/(fs*fs)
  return (obj, obj/dobj )
end

"""
    ssa_simple_newton!(A,grad,PQ,ssa_eps=nothing) -> ssa

Equivalent to `ssa_simple!( ... )`, but uses the Netwon algorightm to compute the SSA.
It might be slightly more efficient. 
"""
function ssa_simple_newton!(A::Matrix{R},
    grad::Union{Nothing,Matrix{R}},
    PQ::SSAAlloc{R,C}, ssa_eps::Union{Nothing,R}=nothing) where {R,C}
  _ssa_eps = something(ssa_eps, default_eps_ssa(A))
  PQ_init!(PQ,A)
  _sa = spectral_abscissa(PQ)
  _start = _sa + 0.1abs(_sa)
  objfun(s)=ssa_simple_obj_newton(s,PQ,_ssa_eps,_sa)
  s_star::Float64 = find_zero(objfun, _start ,Roots.Newton();maxevals=1_000)
  isnothing(grad) && return s_star
  get_P!(s_star,PQ)
  get_Q!(s_star,PQ)
  cc=inv(trace_of_product(PQ.Q,PQ.P))
  mul!(grad,PQ.Q,PQ.P,cc,0.0)
  return s_star
end

#=

  g(s)=__g(s,PQ,stuff, _ssa_eps,_sa)
  dg(s) = __dg(s,PQ,stuff,_sa)
  s_star = find_zero( (g,dg),_start, Roots.Newton() )
  # no grad, end here
  isnothing(grad) && return s_star
  # if grad...
  P_star = _get_P!(s_star,PQ)
  Q_star= _get_Q!(s_star,PQ)
  cc=inv(_trace_prod(Q_star,P_star))
  mul!(grad,Q_star,P_star,cc,0.0)
  return s_star
end

function ssa_simple_newton(A::AbstractArray{R} , ssa_eps::R=zero(R)) where R
    PQ=SSAAlloc(size(A,1))
    return ssa_simple_newton(A,nothing,PQ , ssa_eps)
end

=#

end # module
