var documenterSearchIndex = {"docs":
[{"location":"#Smoothed-Spectral-Abscissa-(SSA)-1","page":"Home","title":"Smoothed Spectral Abscissa (SSA)","text":"","category":"section"},{"location":"#","page":"Home","title":"Home","text":"This package computes the smoothed spectral abscissa (SSA) of square matrices, and the associated gradient, as described in:","category":"page"},{"location":"#","page":"Home","title":"Home","text":"The Smoothed Spectral Abscissa for Robust Stability Optimization , J Vanbiervliet et al , 2009. DOI: 10.1137/070704034","category":"page"},{"location":"#","page":"Home","title":"Home","text":"The SSA is a smooth upper bound to the spectral abscissa of a matrix, that is, the highest real part of the eigenvalues.","category":"page"},{"location":"#","page":"Home","title":"Home","text":"The current version implements only the \"simplified\" version of SSA, i.e. the one where input and output transformations are identity matrices.","category":"page"},{"location":"#Definition-of-SSA-1","page":"Home","title":"Definition of SSA","text":"","category":"section"},{"location":"#","page":"Home","title":"Home","text":"Consider a square real matrix  A, that regulates a dynamical system as follows: textd mathbfxtextdt = Amathbfx. The stability of the dynamical system can be assessed using the following quantity:","category":"page"},{"location":"#","page":"Home","title":"Home","text":"f(As) = int_0^infty textdt left expleft( left(A-sIright)t right)right^2","category":"page"},{"location":"#","page":"Home","title":"Home","text":"where left M right^2 = texttraceleft(M M^top right).","category":"page"},{"location":"#","page":"Home","title":"Home","text":"For a given varepsilon, the SSA can be denoted as tildealpha_varepsilon(A). And satisfies the following  equality:","category":"page"},{"location":"#","page":"Home","title":"Home","text":"f(Atildealpha_varepsilon(A)) = frac1varepsilon","category":"page"},{"location":"#","page":"Home","title":"Home","text":"It is an upper bound to the spectral abscissa of A. Therefore a reduction of SSA to values  0 guarantees dynamical stability.","category":"page"},{"location":"#Interface-1","page":"Home","title":"Interface","text":"","category":"section"},{"location":"#","page":"Home","title":"Home","text":"This module does not export functions in the global scope. It is therefore convenient to shorten the module name as follows:","category":"page"},{"location":"#","page":"Home","title":"Home","text":"using SmoothedSpectralAbscissa ; const SSA=SmoothedSpectralAbscissa","category":"page"},{"location":"#","page":"Home","title":"Home","text":"Therefore one can type SSA.foo in place of SmoothedSpectralAbscissa.foo.","category":"page"},{"location":"#","page":"Home","title":"Home","text":"The functions below compute the SSA (and its gradient) for a matrix A.","category":"page"},{"location":"#","page":"Home","title":"Home","text":"SSA.ssa_simple","category":"page"},{"location":"#SmoothedSpectralAbscissa.ssa_simple","page":"Home","title":"SmoothedSpectralAbscissa.ssa_simple","text":"ssa_simple(A,ssa_eps=nothing) -> ssa\n\nComputes the smoothed spectral abscissa (SSA) of matrix A (with identity input-output weighting matrices). If ssa_eps is not specified, the default value is set by default_eps_ssa(A)\n\nArguments\n\nA::Matrix{<:Real}: A square matrix\nssa_eps::Union{Nothing,R}=nothing The varepsilon parameter associated with the SSA computation. If not specified, it is set by default_eps_ssa(A)\n\nReturns\n\nssa<:Real: this is the tildealpha_varepsilon(A)\n\n\n\n\n\n","category":"function"},{"location":"#","page":"Home","title":"Home","text":"SSA.ssa_simple_withgradient","category":"page"},{"location":"#SmoothedSpectralAbscissa.ssa_simple_withgradient","page":"Home","title":"SmoothedSpectralAbscissa.ssa_simple_withgradient","text":"  ssa_simple_withgradient(A,ssa_eps) -> (ssa,gradient_matrix)\n\nComputes the Smoothed Spectral Abscissa (SSA) of matrix A and also its gradient with respect to each element of A. If ssa_eps is not specified, the default value is set by default_eps_ssa(A)\n\nIn case you need to call the SSA iteratively (e.g. gradient-based optimization), please consider pre-allocating memory and using ssa_simple!(...).\n\n\n\n\n\n","category":"function"},{"location":"#Advanced-Interface-1","page":"Home","title":"Advanced Interface","text":"","category":"section"},{"location":"#","page":"Home","title":"Home","text":"When the SSA is used as optimization objective, it is convenient to use the advanced interface to avoid memory reallocation. The memory is pre-allocated in an object of type SSA.SSAlloc, which can then be used to call the SSA.ssa_simple(...) function multiple times.","category":"page"},{"location":"#","page":"Home","title":"Home","text":"SSA.SSAAlloc","category":"page"},{"location":"#SmoothedSpectralAbscissa.SSAAlloc","page":"Home","title":"SmoothedSpectralAbscissa.SSAAlloc","text":"    SSAAlloc(n::Integer) -> SSAAlloc\n\nMemory allocation to avoid re-allocating space for repeated SSA computations.\n\nArguments\n\nn::Integer: the size of the matrix (or matrices) on which the SSA will be computed\n\n\n\n\n\n    SSAAlloc(A::Matrix) -> SSAAlloc\n\nMemory allocation to avoid re-allocating space for repeated SSA computations. This corresponds to SSAAlloc(size(A,1))\n\n\n\n\n\n","category":"type"},{"location":"#","page":"Home","title":"Home","text":"Once the space is allocated, the SSA and its gradient can be computed by the following functions","category":"page"},{"location":"#","page":"Home","title":"Home","text":"SSA.ssa_simple!","category":"page"},{"location":"#SmoothedSpectralAbscissa.ssa_simple!","page":"Home","title":"SmoothedSpectralAbscissa.ssa_simple!","text":"ssa_simple!(A,grad,PQ,ssa_eps=nothing) -> ssa\n\nComputes the smoothed spectral abscissa (SSA) of matrix A (with identity input-output weighting matrices) and its gradient (optionally). For efficiency, A is reallocated in the computation. If ssa_eps is not specified, the default value is set by default_eps_ssa(A)\n\nArguments\n\nA::Matrix{<:Real}: A square matrix\ngrad::Union{Nothing,Matrix{<:Real}}: nothing  skips the gradient computation. To compute the gradient, provide a matrix grad of the same size as A. The matrix will be rewritten in-place with entries fracpartial tildealpha_varepsilon(A)partial A_ij\nPQ::SSAAlloc : Pre-allocated matrices for iterative computations. Can be generated by SSAAlloc(A)\nssa_eps::Union{Nothing,R}=nothing The varepsilon parameter associated with the SSA computation. If not specified, it is set by default_eps_ssa(A)\n\nReturns\n\nssa<:Real: this is tildealpha_varepsilon(A)\n\n\n\n\n\n","category":"function"},{"location":"#","page":"Home","title":"Home","text":"SSA.ssa_simple_newton!","category":"page"},{"location":"#SmoothedSpectralAbscissa.ssa_simple_newton!","page":"Home","title":"SmoothedSpectralAbscissa.ssa_simple_newton!","text":"ssa_simple_newton!(A,grad,PQ,ssa_eps=nothing) -> ssa\n\nEquivalent to ssa_simple!( ... ), but uses the Netwon algorightm to compute the SSA. It might be slightly more efficient. \n\n\n\n\n\n","category":"function"},{"location":"#Index-1","page":"Home","title":"Index","text":"","category":"section"},{"location":"#","page":"Home","title":"Home","text":"","category":"page"}]
}
