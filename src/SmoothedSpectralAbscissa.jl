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


"""
    Workspace(n::Integer)
    Workspace(A::Matrix{<:Real})

Allocate reusable storage for SSA computations on nonempty square matrices.
The twelve working matrices use `Float64`; cached eigenvalues use `ComplexF32`.
`Workspace(A)` uses the size of `A`, irrespective of its element type, without
copying it. Storage is uninitialized until `ssa` or `PQ_init!` is called.

The type parameters describe the concrete matrix and vector container types:
`Workspace{Matrix{Float64},Vector{ComplexF32}}` for these constructors.
A workspace can be reused for different matrices of the same size. It is
mutated during computation and must not be shared by concurrent calls.
"""
struct Workspace{M,V}
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
    A_eigvals::V
end

function Workspace(n::Integer)
    if n <= 0
        throw(ArgumentError("Workspace size must be positive."))
    end
    return Workspace(
        map(_ -> Matrix{Float64}(undef, n, n), 1:12)...,
        Vector{ComplexF32}(undef, n))
end

function Workspace(A::Matrix{<:Real})
    if size(A, 1) != size(A, 2)
        throw(DimensionMismatch("SSA requires a square matrix."))
    end
    return Workspace(size(A, 1))
end

"""
    PQ_init!(PQ::Workspace{Matrix{R},V}, A::Matrix{R}) -> nothing

Initialize `PQ` with Schur decompositions of `A` and its transpose, the
eigenvalues of `A`, and the transformed identity weighting matrices.
Call once for each new `A`, before starting the root-finding procedure.
`ssa` calls this function internally.

# Arguments
- `PQ::Workspace{Matrix{R},V}`: Working storage with matrices of the same size and element
  type as `A`. Its cached decompositions and eigenvalues are overwritten.
- `A::Matrix{R}`: Real square matrix to decompose. `A` is unchanged.

# Returns
- `nothing`: The precomputed quantities are stored in `PQ`.
"""
function PQ_init!(PQ::Workspace{Matrix{R},V},A::Matrix{R}) where {R,V}
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

@inline function spectral_abscissa(PQ::Workspace)
    return convert(eltype(PQ.R), maximum(real.(PQ.A_eigvals)))
end
function get_P_or_Q!(PorQ::M,s::Real,Z::M,R::M,D::M,
        YZt_alloc::M) where M<:Matrix{<:Real}
    subtract_diagonal!(R,s)
    Y, scale = LAPACK.trsyl!('N','T', R, R, D)
    mul!(YZt_alloc,Y,transpose(Z))
    mul!(PorQ,Z,YZt_alloc,inv(scale),0.0)
    return nothing
end
function get_P!(s::T, PQ::Workspace{Matrix{T},V}) where {T,V}
    R=copyto!(PQ.R_alloc,PQ.R)
    D=copyto!(PQ.D_alloc,PQ.D)
    get_P_or_Q!(PQ.P,s,PQ.Z,R,D,PQ.YZt_alloc)
    return nothing
end
function get_Q!(s::T, PQ::Workspace{Matrix{T},V}) where {T,V}
    R=copyto!(PQ.R_alloc,PQ.Rt)
    D=copyto!(PQ.D_alloc,PQ.Dt)
    get_P_or_Q!(PQ.Q,s,PQ.Zt,R,D,PQ.YZt_alloc)
    return nothing
end

function ssa_simple_obj(s::R,PQ::Workspace{Matrix{R},V},ssa_eps::R,sa::R) where {R,V}
    get_P!( max(sa+eps(10.0),s) ,PQ) # SSA is not defined below SA !
    return inv(tr(PQ.P)) - ssa_eps
end

# Find zero with Newton!
function ssa_simple_obj_newton(s::R,PQ::Workspace{Matrix{R},V},ssa_eps::R,sa::R) where {R,V}
  _s = max(sa+eps(10.0),s)
  get_P!(_s,PQ) # SSA is not defined below SA !
  get_Q!(_s,PQ)
  fs = tr(PQ.P)
  mdfs = 2.0*trace_of_product(PQ.P,PQ.Q)
  obj = inv(fs) - ssa_eps
  dobj = mdfs/(fs*fs)
  return (obj,obj/dobj)
end


function validate_ssa_inputs(A::Matrix{Float64}, grad, workspace::Workspace,
        ssa_eps::Float64)
    if size(A, 1) != size(A, 2)
        throw(DimensionMismatch("SSA requires a square matrix."))
    end
    if isempty(A)
        throw(ArgumentError("SSA requires a nonempty matrix."))
    end
    if !isfinite(ssa_eps) || ssa_eps <= 0
        throw(ArgumentError("The smoothing parameter must be finite and positive."))
    end
    if !isnothing(grad)
        if size(grad) != size(A)
            throw(DimensionMismatch("Gradient dimensions must match A."))
        end
        if Base.mightalias(A, grad)
            throw(ArgumentError("The gradient must not alias A."))
        end
    end
    for name in fieldnames(typeof(workspace))
        buffer = getfield(workspace, name)
        expected_size = size(A)
        if name === :A_eigvals
            expected_size = (size(A, 1),)
        end
        if size(buffer) != expected_size
            throw(DimensionMismatch("Workspace buffer $name has incompatible dimensions."))
        end
        if Base.mightalias(A, buffer)
            throw(ArgumentError("Workspace buffers must not alias A."))
        end
        if !isnothing(grad)
            if Base.mightalias(grad, buffer)
                throw(ArgumentError("The gradient must not alias workspace buffers."))
            end
        end
    end
    return nothing
end

"""
    ssa(A, ssa_eps=nothing; workspace=nothing, grad=nothing,
        optim_method=SmoothedSpectralAbscissa.OptimOrder2,
        input_output_weighting=LinearAlgebra.I) -> ssa_value

Compute the smoothed spectral abscissa (SSA) of a nonempty square
`Matrix{Float64}`. `A` is unchanged. Despite the absence of a `!` suffix,
`workspace` and an optional `grad::Matrix{Float64}` are overwritten.

For identity input/output weighting, the SSA is the shift `s` above the spectral
abscissa of `A` for which `inv(tr(P)) == ssa_eps`, where `P` solves
`(A - s*I)*P + P*(A - s*I)' + I = 0`.

`ssa_eps` must be a finite positive `Float64`; `nothing` uses
`default_eps_ssa(A)`. Pass `workspace=Workspace(A)` to reuse working storage;
when omitted or `nothing`, storage is allocated internally. Workspace dimensions
must match `A` and its matrix buffers must use `Float64`.

Pass `grad=similar(A)` to store the gradient while returning only the SSA value.
Entry `grad[i, j]` is the derivative with respect to `A[i, j]`. The gradient
must match the dimensions of `A`. Input, gradient, and workspace storage must
not alias each other.

`optim_method` accepts the types `OptimOrder2` (default) and `OptimNewton`,
not instances. Only `input_output_weighting=LinearAlgebra.I` is implemented;
explicit identity matrices and other uniform scaling operators raise an error.
"""
function ssa(A::Matrix{Float64}, ssa_eps::Union{Nothing,Float64}=nothing;
        workspace::Union{Nothing,Workspace{Matrix{Float64}}}=nothing,
        grad::Union{Nothing,Matrix{Float64}}=nothing,
        optim_method::Type=OptimOrder2,
        input_output_weighting::Union{UniformScaling,AbstractMatrix}=I)
    if input_output_weighting !== I
        error("Input-output weighting other than LinearAlgebra.I is not implemented.")
    end
    if isnothing(workspace)
        workspace = Workspace(A)
    end
    epsilon = something(ssa_eps, default_eps_ssa(A))
    validate_ssa_inputs(A, grad, workspace, epsilon)
    return _ssa!(A, grad, workspace, epsilon, optim_method)
end

# Order2() to find roots
function _ssa!(A::Matrix{R},grad::Union{Nothing,Matrix{R}}, alloc::Workspace{Matrix{R},V},
     ssa_eps::R , optim_method::Type{OptimOrder2}) where {R,V}
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
function _ssa!(A::Matrix{R},grad::Union{Nothing,Matrix{R}}, alloc::Workspace{Matrix{R},V},
     ssa_eps::R , optim_method::Type{OptimNewton}) where {R,V}
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
    ssa_withgradient(A, ssa_eps=nothing; workspace=nothing,
        optim_method=SmoothedSpectralAbscissa.OptimOrder2,
        input_output_weighting=LinearAlgebra.I) -> (ssa_value, gradient)

Compute the SSA and allocate its gradient matrix. Accepts the same input,
smoothing parameter, workspace, solver, and weighting options as [`ssa`](@ref).
`A` is unchanged; a supplied workspace is overwritten and reused. A workspace
is allocated internally when omitted or `nothing`.

For repeated computations, preallocate both a workspace and gradient matrix,
then call `ssa(A, ssa_eps; workspace=workspace, grad=gradient)`.
"""
function ssa_withgradient(A::Matrix{Float64},
        ssa_eps::Union{Nothing,Float64}=nothing;
        workspace::Union{Nothing,Workspace{Matrix{Float64}}}=nothing,
        optim_method::Type=OptimOrder2,
        input_output_weighting::Union{UniformScaling,AbstractMatrix}=I)
    gradient = similar(A)
    value = ssa(A, ssa_eps; workspace=workspace, grad=gradient,
        optim_method=optim_method, input_output_weighting=input_output_weighting)
    return value, gradient
end

end # module
