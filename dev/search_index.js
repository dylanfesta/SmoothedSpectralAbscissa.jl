var documenterSearchIndex = {"docs":
[{"location":"#SmoothedSpectralAbscissa.jl-1","page":"Home","title":"SmoothedSpectralAbscissa.jl","text":"","category":"section"},{"location":"#","page":"Home","title":"Home","text":"The Purpose of this package is computing the Smoothed Spectral Abscissa of square matrices, and the associated gradient. The computation can be optimized for recursion, so that each new calculation does not reallocate memory.","category":"page"},{"location":"#","page":"Home","title":"Home","text":"I implemented only the \"simplified\" version that does not have projections.","category":"page"},{"location":"#","page":"Home","title":"Home","text":"The algorithm is described in the following paper:","category":"page"},{"location":"#","page":"Home","title":"Home","text":"The Smoothed Spectral Abscissa for Robust Stability Optimization , J Vanbiervliet et al , 2009. DOI: 10.1137/070704034","category":"page"},{"location":"#TO-DOs-1","page":"Home","title":"TO-DOs","text":"","category":"section"},{"location":"#","page":"Home","title":"Home","text":"Add references\ndescribe the algorithm   ","category":"page"},{"location":"#","page":"Home","title":"Home","text":"","category":"page"},{"location":"#","page":"Home","title":"Home","text":"","category":"page"},{"location":"#","page":"Home","title":"Home","text":"Modules = [SmoothedSpectralAbscissa]","category":"page"},{"location":"#SmoothedSpectralAbscissa.SSAAlloc-Tuple{Integer}","page":"Home","title":"SmoothedSpectralAbscissa.SSAAlloc","text":"    SSAAlloc{M<:AbstractMatrix,V<:AbstractVector}(n::Integer) -> SSAAlloc\n\nAllocation structure that makes SSA computation mode efficient. n is the size of the matrix (or matrices) on which one wants to compute the SSA\n\n\n\n\n\n","category":"method"},{"location":"#SmoothedSpectralAbscissa.PQ_init!-Tuple{AbstractArray{T,2} where T,SmoothedSpectralAbscissa.SSAAlloc}","page":"Home","title":"SmoothedSpectralAbscissa.PQ_init!","text":"    PQ_init!(A::AbstractMatrix,PQ::SSAAlloc) -> nothing\n\nPrepares SSAAlloc with the matrix on which one wants to compute the SSA. Schur decompositions are performed here.\n\n\n\n\n\n","category":"method"},{"location":"#SmoothedSpectralAbscissa.default_eps_ssa-Tuple{AbstractArray{T,2} where T}","page":"Home","title":"SmoothedSpectralAbscissa.default_eps_ssa","text":"    default_eps_ssa(A::AbstractMatrix{Float64}) -> eps_ssa\n\ndefault value for the varepsilon used to compute the SSA\n\n\n\n\n\n","category":"method"},{"location":"#SmoothedSpectralAbscissa.ssa_simple","page":"Home","title":"SmoothedSpectralAbscissa.ssa_simple","text":"    ssa_simple(A::AbstractMatrix , ssa_eps::Union{Nothing,Float64}=nothing)\n        -> ssa::Float64\n\nComputes the Smoothed Spectral Abscissa (SSA) of matrix A in its simplest formulation If ssa_eps is not specified, a default value is used.\n\nIn case you need to call the SSA iteratively, consider the version that acts on pre-allocated memory.\n\n\n\n\n\n","category":"function"},{"location":"#SmoothedSpectralAbscissa.ssa_simple!","page":"Home","title":"SmoothedSpectralAbscissa.ssa_simple!","text":"    ssa_simple!(A::AbstractMatrix, with_gradient::Bool,PQ::SSAAlloc ,\n        ssa_eps::Union{Nothing,Float64}=nothing) -> ssa::Float64\n\nComputes the Smoothed Spectral Abscissa (SSA) of matrix A in its simplest formulation. For efficiency, A is reallocated in the computation. If ssa_eps is not specified, a default value is used.\n\n\n\n\n\n","category":"function"},{"location":"#SmoothedSpectralAbscissa.ssa_simple_newton","page":"Home","title":"SmoothedSpectralAbscissa.ssa_simple_newton","text":"This uses the Newton method and the analytic derivative , seems as fast as the Order16 method\n\n\n\n\n\n","category":"function"},{"location":"#SmoothedSpectralAbscissa.ssa_simple_withgradient","page":"Home","title":"SmoothedSpectralAbscissa.ssa_simple_withgradient","text":"  ssa_simple_withgradient(A::AbstractMatrix , ssa_eps::Union{Nothing,Float64}=nothing)\n        -> (ssa,gradient_matrix)::Tuple{Float64,Matrix}\n\nComputes the Smoothed Spectral Abscissa (SSA) of matrix A in its simplest formulation, and also its gradient, for each matrix element. If ssa_eps is not specified, a default value is used.\n\nIn case you need to call the SSA iteratively (e.g. gradient-based optimization),  consider the version that acts on pre-allocated memory.\n\n\n\n\n\n","category":"function"}]
}
