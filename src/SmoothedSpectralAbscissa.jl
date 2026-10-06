module SmoothedSpectralAbscissa
import LinearAlgebra: I, LAPACK, UniformScaling, eigvals, mul!, schur, schur!,
    tr, transpose!
import Roots: Newton, Order2, find_zero


abstract type OptimMethod end
struct OptimOrder2 end
struct OptimNewton end


@inline function spectral_abscissa(A::Matrix{R}) where R
   convert(R,maximum(real.(eigvals(A))))
end

# trace of matrix product
@inline function trace_of_product(A::Matrix{R},B::Matrix{R}) where R
    n=size(A,1)
    acc=zero(R)
    for i in 1:n, j in 1:n
        @inbounds acc += A[i,j]*B[j,i]
    end
    return acc
end
# remove value from the diagonal
# equivalent to copy!(A,A-s*I)
@inline function subtract_diagonal!(A::Matrix{R},s::R) where R
    n=size(A,1)
    @simd for i in 1:n
        @inbounds A[i,i] -= s
    end
    return nothing
end

"""
    default_eps_ssa(A::Matrix{<:Real}) -> ssa_eps

Return the default smoothing parameter for the smoothed spectral abscissa (SSA).
It scales inversely with the number of rows of `A` and equals `0.01` for a
``150 \\times 150`` matrix.

# Arguments
- `A::Matrix{<:Real}`: Matrix whose number of rows sets the smoothing parameter.

# Returns
- `ssa_eps::Float64`: Default smoothing parameter, `0.01 * 150.0 / size(A, 1)`.
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
    SSAAlloc(n::Integer) -> alloc

Allocate reusable working storage for SSA computations on `n`-by-`n` matrices.
The working matrices use `Float64`; cached eigenvalues use `ComplexF32`.

# Arguments
- `n::Integer`: Number of rows and columns of the matrices used in SSA
  computations.

# Returns
- `alloc::SSAAlloc{Float64,ComplexF32}`: Uninitialized working storage. `ssa!`
  initializes the storage for the supplied matrix before computing the SSA.
"""
function SSAAlloc(n::Integer)
  return SSAAlloc(
   map( _ -> Matrix{Float64}(undef,n,n),1:12)... , Vector{ComplexF32}(undef,n))
end

"""
    SSAAlloc(A::Matrix{<:Real}) -> alloc

Allocate reusable working storage with `SSAAlloc(size(A, 1))`.
The element type of `A` does not change the storage types, and this constructor
neither copies `A` nor initializes its Schur decompositions.

# Arguments
- `A::Matrix{<:Real}`: Square matrix whose size determines the working storage.

# Returns
- `alloc::SSAAlloc{Float64,ComplexF32}`: Uninitialized working storage with
  `size(A, 1)`-by-`size(A, 1)` working matrices.
"""
function SSAAlloc(A::Matrix{<:Real})
  return SSAAlloc(size(A,1))
end

"""
    PQ_init!(PQ::SSAAlloc{R,C}, A::Matrix{R}) -> nothing

Initialize `PQ` with Schur decompositions of `A` and its transpose, the
eigenvalues of `A`, and the transformed identity weighting matrices.
Call once for each new `A`, before starting the root-finding procedure.
`ssa!` calls this function internally.

# Arguments
- `PQ::SSAAlloc{R,C}`: Working storage with matrices of the same size and element
  type as `A`. Its cached decompositions and eigenvalues are overwritten.
- `A::Matrix{R}`: Real square matrix to decompose. `A` is unchanged.

# Returns
- `nothing`: The precomputed quantities are stored in `PQ`.
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
    get_P!( max(sa+eps(10.0),s) ,PQ) # SSA is not defined below SA !
    return inv(tr(PQ.P)) - ssa_eps
end

# Find zero with Newton!
function ssa_simple_obj_newton(s::R,PQ::SSAAlloc{R,C},ssa_eps::R,sa::R) where {R,C}
  _s = max(sa+eps(10.0),s)
  get_P!(_s,PQ) # SSA is not defined below SA !
  get_Q!(_s,PQ)
  fs = tr(PQ.P)
  mdfs = 2.0*trace_of_product(PQ.P,PQ.Q)
  obj = inv(fs) - ssa_eps
  dobj = mdfs/(fs*fs)
  return (obj,obj/dobj)
end


"""
    ssa!(A, grad, alloc, ssa_eps=nothing;
         optim_method=SmoothedSpectralAbscissa.OptimOrder2,
         input_output_weighting=LinearAlgebra.I) -> ssa_value

Compute the smoothed spectral abscissa (SSA) of `A` and optionally its gradient,
using reusable working storage. `A` is unchanged; `alloc` and, when supplied,
`grad` are overwritten.

For identity input/output weighting, the SSA is the shift `s` above the spectral
abscissa of `A` for which `inv(tr(P)) == ssa_eps`, where `P` solves
`(A - s*I)*P + P*(A - s*I)' + I = 0`.

# Arguments
- `A::Matrix{R}`: Real square matrix. With the current `SSAAlloc` constructors,
  `R` must be `Float64`.
- `grad::Union{Nothing,Matrix{R}}`: Pass `nothing` to skip the gradient computation,
  or a matrix of the same size and element type as `A` to store the gradient.
  Entry `grad[i, j]` is the derivative of the SSA with respect to `A[i, j]`.
- `alloc::SSAAlloc{R,C}`: Reusable working storage with matrices of the same size
  and element type as `A`, created by `SSAAlloc(A)` or `SSAAlloc(size(A, 1))`.
- `ssa_eps::Union{Nothing,R}=nothing`: Positive smoothing parameter
  ``\\varepsilon``. Pass `nothing` to use `default_eps_ssa(A)`.

# Keyword Arguments
- `optim_method::Type=SmoothedSpectralAbscissa.OptimOrder2`: Root-finding method.
  The supported types are `OptimOrder2` and `OptimNewton` from this module.
  Pass the type itself rather than an instance.
- `input_output_weighting::Union{UniformScaling,AbstractMatrix}=LinearAlgebra.I`:
  Input/output weighting. Only `LinearAlgebra.I` is implemented. Matrices and
  other uniform scaling operators throw an error indicating that the weighting
  is not implemented, including explicit identity matrices.

# Returns
- `ssa_value::R`: Smoothed spectral abscissa ``\\tilde{\\alpha}_\\varepsilon(A)``.
  The gradient, when requested, is stored in `grad`.
"""
function ssa!(A::Matrix{R},grad::Union{Nothing,Matrix{R}}, alloc::SSAAlloc{R,C},
    ssa_eps::Union{Nothing,R}=nothing ; optim_method::Type=OptimOrder2,
    input_output_weighting::Union{UniformScaling,AbstractMatrix}=I) where {R,C}
  if input_output_weighting !== I
    error("Input-output weighting other than LinearAlgebra.I is not implemented.")
  end
  _ssa_eps = something(ssa_eps, default_eps_ssa(A))
  return ssa!(A,grad,alloc,_ssa_eps,optim_method)
end
# Order2() to find roots
function ssa!(A::Matrix{R},grad::Union{Nothing,Matrix{R}}, alloc::SSAAlloc{R,C},
     ssa_eps::R , optim_method::Type{OptimOrder2}) where {R,C}
  PQ_init!(alloc,A)
  _sa = spectral_abscissa(alloc)
  _start = _sa + 0.5*ssa_eps
  objfun(s)=ssa_simple_obj(s,alloc,ssa_eps,_sa)
  _method = Order2()
  s_star::R = find_zero(objfun, _start,_method;maxevals=1_000)
  if isnothing(grad)
    return s_star
  end
  # compute gradient
  get_P!(s_star,alloc)
  get_Q!(s_star,alloc)
  cc=inv(trace_of_product(alloc.Q,alloc.P))
  mul!(grad,alloc.Q,alloc.P,cc,0.0)
  return s_star
end

# newton to find roots
function ssa!(A::Matrix{R},grad::Union{Nothing,Matrix{R}}, alloc::SSAAlloc{R,C},
     ssa_eps::R , optim_method::Type{OptimNewton}) where {R,C}
  PQ_init!(alloc,A)
  _sa = spectral_abscissa(alloc)
  _start = _sa + 0.5*ssa_eps # _start = _sa + 0.1abs(_sa)
  objfun(s)=ssa_simple_obj_newton(s,alloc,ssa_eps,_sa)
  _method = Newton()
  s_star::R = find_zero(objfun, _start,_method;maxevals=1_000)
  if isnothing(grad)
    return s_star
  end
  # compute gradient
  get_P!(s_star,alloc)
  get_Q!(s_star,alloc)
  cc=inv(trace_of_product(alloc.Q,alloc.P))
  mul!(grad,alloc.Q,alloc.P,cc,0.0)
  return s_star
end


"""
    ssa_simple!(A, grad, PQ, ssa_eps=nothing) -> ssa_value

Legacy wrapper retained for testing. Equivalent to `ssa!(A, grad, PQ, ssa_eps)`
with the default root-finding method and identity input/output weighting.
Use `ssa!` for new code.

# Arguments
- `A::Matrix{R}`: Real square matrix. With the current `SSAAlloc` constructors,
  `R` must be `Float64`.
- `grad::Union{Nothing,Matrix{R}}`: Pass `nothing` to skip the gradient computation,
  or a matrix of the same size and element type as `A` to store the gradient.
- `PQ::SSAAlloc{R,C}`: Reusable working storage with matrices of the same size
  and element type as `A`. Overwritten during the computation.
- `ssa_eps::Union{Nothing,R}=nothing`: Positive smoothing parameter
  ``\\varepsilon``. Pass `nothing` to use `default_eps_ssa(A)`.

# Returns
- `ssa_value::R`: Smoothed spectral abscissa ``\\tilde{\\alpha}_\\varepsilon(A)``.
  The gradient, when requested, is stored in `grad`; `A` is unchanged.
"""
function ssa_simple!(A::Matrix{R},grad::Union{Nothing,Matrix{R}},
        PQ::SSAAlloc{R,C},ssa_eps::Union{Nothing,R}=nothing) where {R,C}
 return ssa!(A,grad,PQ,ssa_eps)
end


# short versions that also allocates the memory
"""
    ssa(A, ssa_eps=nothing,
        input_output_weighting=LinearAlgebra.I) -> ssa_value

Compute the smoothed spectral abscissa (SSA) of `A`, allocating working storage
internally. `A` is unchanged. See `ssa!` for the defining equation.

# Arguments
- `A::Matrix{Float64}`: Real square matrix. The current working-storage
  constructors require `Float64` inputs for this computation.
- `ssa_eps::Union{Nothing,Float64}=nothing`: Positive smoothing parameter
  ``\\varepsilon``. Pass `nothing` to use `default_eps_ssa(A)`.
- `input_output_weighting::Union{UniformScaling,AbstractMatrix}=LinearAlgebra.I`:
  Optional third positional argument specifying input/output weighting. Only
  `LinearAlgebra.I` is implemented. Matrices and other uniform scaling operators
  throw an error indicating that the weighting is not implemented, including
  explicit identity matrices.

# Returns
- `ssa_value::Float64`: Smoothed spectral abscissa
  ``\\tilde{\\alpha}_\\varepsilon(A)``.
"""
function ssa(A::Matrix{R},ssa_eps::Union{Nothing,R}=nothing,
    input_output_weighting::Union{UniformScaling,AbstractMatrix}=I) where R
  alloc=SSAAlloc(size(A,1))
  _epsssa = something(ssa_eps, default_eps_ssa(A))
  return ssa!(copy(A),nothing,alloc, _epsssa;input_output_weighting=input_output_weighting)
end

"""
    ssa_withgradient(A, ssa_eps=nothing,
                     input_output_weighting=LinearAlgebra.I)
        -> (ssa_value, gradmat)

Compute the smoothed spectral abscissa (SSA) of `A` and its gradient with respect
to each entry of `A`, allocating working storage and the gradient matrix
internally. `A` is unchanged. See `ssa!` for the defining equation.

For repeated computations, preallocate working storage and a gradient matrix,
then use `ssa!`.

# Arguments
- `A::Matrix{Float64}`: Real square matrix. The current working-storage
  constructors require `Float64` inputs for this computation.
- `ssa_eps::Union{Nothing,Float64}=nothing`: Positive smoothing parameter
  ``\\varepsilon``. Pass `nothing` to use `default_eps_ssa(A)`.
- `input_output_weighting::Union{UniformScaling,AbstractMatrix}=LinearAlgebra.I`:
  Optional third positional argument specifying input/output weighting. Only
  `LinearAlgebra.I` is implemented. Matrices and other uniform scaling operators
  throw an error indicating that the weighting is not implemented, including
  explicit identity matrices.

# Returns
A tuple `(ssa_value, gradmat)` containing:

- `ssa_value::Float64`: Smoothed spectral abscissa
  ``\\tilde{\\alpha}_\\varepsilon(A)``.
- `gradmat::Matrix{Float64}`: Gradient matrix of the same size as `A`.
  Entry `gradmat[i, j]` is the derivative of the SSA with respect to `A[i, j]`.
"""
function ssa_withgradient(A::Matrix{R},ssa_eps::Union{Nothing,R}=nothing,
        input_output_weighting::Union{UniformScaling,AbstractMatrix}=I) where R
    alloc=SSAAlloc(size(A,1))
    gradmat=similar(A)
    epsssa = something(ssa_eps, default_eps_ssa(A))
    _ssa = ssa!(copy(A),gradmat,alloc, epsssa;input_output_weighting=input_output_weighting)
    return _ssa,gradmat
end


end # module
