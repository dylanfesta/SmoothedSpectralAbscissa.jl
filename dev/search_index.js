var documenterSearchIndex = {"docs":
[{"location":"02_dynamics/","page":"Faster convergence for linear dynamcs","title":"Faster convergence for linear dynamcs","text":"EditURL = \"https://github.com/dylanfesta/SmoothedSpectralAbscissa.jl/blob/master/examples/02_dynamics.jl\"","category":"page"},{"location":"02_dynamics/#Faster-convergence-for-linear-dynamcs","page":"Faster convergence for linear dynamcs","title":"Faster convergence for linear dynamcs","text":"","category":"section"},{"location":"02_dynamics/","page":"Faster convergence for linear dynamcs","title":"Faster convergence for linear dynamcs","text":"In this example, I show how the SSA can be used as objective to guarantee convergence or faster convergence for a linear dynamical system.","category":"page"},{"location":"02_dynamics/#Optimization-of-all-matrix-elements","page":"Faster convergence for linear dynamcs","title":"Optimization of all matrix elements","text":"","category":"section"},{"location":"02_dynamics/#Initialization","page":"Faster convergence for linear dynamcs","title":"Initialization","text":"","category":"section"},{"location":"02_dynamics/","page":"Faster convergence for linear dynamcs","title":"Faster convergence for linear dynamcs","text":"using Plots,NamedColors\nusing LinearAlgebra\nusing Random\nusing SmoothedSpectralAbscissa ; const SSA=SmoothedSpectralAbscissa\nRandom.seed!(0);\nnothing #hide","category":"page"},{"location":"02_dynamics/#Linear-Dynamics","page":"Faster convergence for linear dynamcs","title":"Linear Dynamics","text":"","category":"section"},{"location":"02_dynamics/","page":"Faster convergence for linear dynamcs","title":"Faster convergence for linear dynamcs","text":"Consider continuous linear dynamics, regulated by","category":"page"},{"location":"02_dynamics/","page":"Faster convergence for linear dynamcs","title":"Faster convergence for linear dynamcs","text":"fractextd mathbfxtextdt = A mathbfx","category":"page"},{"location":"02_dynamics/","page":"Faster convergence for linear dynamcs","title":"Faster convergence for linear dynamcs","text":"The soluton is analytic and takes the form:","category":"page"},{"location":"02_dynamics/","page":"Faster convergence for linear dynamcs","title":"Faster convergence for linear dynamcs","text":"mathbfx(t) = exp( At )mathbfx_0","category":"page"},{"location":"02_dynamics/","page":"Faster convergence for linear dynamcs","title":"Faster convergence for linear dynamcs","text":"Where mathbfx_0 are the initial conditions. The function below computes the dynamical evolution of the system.","category":"page"},{"location":"02_dynamics/","page":"Faster convergence for linear dynamcs","title":"Faster convergence for linear dynamcs","text":"function run_linear_dyn(A::Matrix{R},x0::Vector{R},tmax::Real,dt::Real=0.01) where R\n  ts = range(0,tmax;step=dt)\n  ret = Matrix{R}(undef,length(x0),length(ts))\n  for (k,t) in enumerate(ts)\n    ret[:,k]=exp(A.*t)*x0\n  end\n  retnrm = mapslices(norm,ret;dims=1)[:]\n  return ts,ret,retnrm\nend;\nnothing #hide","category":"page"},{"location":"02_dynamics/#The-optimization-of-the-objective-is-done-through-Optim.jl-and-BFGS","page":"Faster convergence for linear dynamcs","title":"The optimization of the objective is done through Optim.jl and BFGS","text":"","category":"section"},{"location":"02_dynamics/","page":"Faster convergence for linear dynamcs","title":"Faster convergence for linear dynamcs","text":"using Optim\n\nfunction objective_and_grad_simple(x::Vector{R},grad::Union{Nothing,Vector{R}},\n    n::Integer,ssa_eps::R,alloc::SSA.SSAAlloc) where R\n  mat=reshape(x,(n,n))\n  gradmat = isnothing(grad) ? nothing : similar(mat)\n  obj = SSA.ssa_simple!(mat,gradmat,alloc,ssa_eps)\n  if !isnothing(grad)\n    for i in eachindex(gradmat)\n      grad[i]=gradmat[i]\n    end\n  end\n  return obj\nend;\nnothing #hide","category":"page"},{"location":"02_dynamics/#Start-with-an-unstable-matrix","page":"Faster convergence for linear dynamcs","title":"Start with an unstable matrix","text":"","category":"section"},{"location":"02_dynamics/","page":"Faster convergence for linear dynamcs","title":"Faster convergence for linear dynamcs","text":"n = 50\nA = randn(n,n) ./ sqrt(n) + 0.2I\nx0=randn(n)\ntimes,_,dyn_norms = run_linear_dyn(A,x0,3.,0.1)\n\nplot(times,dyn_norms; leg=false,linewidth=3,color=:black,xlabel=\"time\",ylabel=\"norm(x(t))\")","category":"page"},{"location":"02_dynamics/","page":"Faster convergence for linear dynamcs","title":"Faster convergence for linear dynamcs","text":"as expected, the norm grows exponentially with time.","category":"page"},{"location":"02_dynamics/#Now-do-the-gradient-based-optimization","page":"Faster convergence for linear dynamcs","title":"Now do the gradient-based optimization","text":"","category":"section"},{"location":"02_dynamics/","page":"Faster convergence for linear dynamcs","title":"Faster convergence for linear dynamcs","text":"const ssa_eps=0.001\nconst alloc = SSA.SSAAlloc(n)\nconst y0 = A[:];\nnothing #hide","category":"page"},{"location":"02_dynamics/","page":"Faster convergence for linear dynamcs","title":"Faster convergence for linear dynamcs","text":"The objective function to be minimized is","category":"page"},{"location":"02_dynamics/","page":"Faster convergence for linear dynamcs","title":"Faster convergence for linear dynamcs","text":"textobj(A) = textSSA(A) + lambda frac12 left A - A_0 right^2","category":"page"},{"location":"02_dynamics/","page":"Faster convergence for linear dynamcs","title":"Faster convergence for linear dynamcs","text":"Where lambda sets the relative weight. We are reducing the SSA while keeping the matrix elements close to their initial value. Adding some form of regularization is always necessary when otpimizing. If not, the SSA would run to -infty.","category":"page"},{"location":"02_dynamics/","page":"Faster convergence for linear dynamcs","title":"Faster convergence for linear dynamcs","text":"function objfun!(F,G,y)\n  λ = 50.0/length(y) # regularizer weight\n  obj=objective_and_grad_simple(y,G,n,ssa_eps,alloc) # SSA and gradient\n  ydiffs = y.-y0\n  obj += 0.5*λ*mapreduce(x->x^2,+,ydiffs) # add the regularizer\n  if !isnothing(G)\n    @. G += λ*ydiffs # add gradient of regularizer\n  end\n  return obj\nend;\nnothing #hide","category":"page"},{"location":"02_dynamics/#Optimize-and-show-the-results","page":"Faster convergence for linear dynamcs","title":"Optimize and show the results","text":"","category":"section"},{"location":"02_dynamics/","page":"Faster convergence for linear dynamcs","title":"Faster convergence for linear dynamcs","text":"opt_out = optimize(Optim.only_fg!(objfun!),A[:],BFGS(),Optim.Options(iterations=50))\ny_opt=Optim.minimizer(opt_out)\nA_opt = reshape(y_opt,(n,n))\ntimes,_,dyn_norms_opt = run_linear_dyn(A_opt,x0,3.,0.1)\nplot(times,dyn_norms_opt;\n  leg=false,linewidth=3,color=:blue,xlabel=\"time\",ylabel=\"norm(x(t)) optimized\")","category":"page"},{"location":"02_dynamics/","page":"Faster convergence for linear dynamcs","title":"Faster convergence for linear dynamcs","text":"The optimized matrix produces stable dynamics, as shown in the plot above.","category":"page"},{"location":"02_dynamics/","page":"Faster convergence for linear dynamcs","title":"Faster convergence for linear dynamcs","text":"We can also take a look at the matrices before and after optimization","category":"page"},{"location":"02_dynamics/","page":"Faster convergence for linear dynamcs","title":"Faster convergence for linear dynamcs","text":"heatmap(hcat(A,fill(NaN,n,10),A_opt);ratio=1,axis=nothing,ticks=nothing,border=:none,\n  colorbar=nothing)","category":"page"},{"location":"02_dynamics/","page":"Faster convergence for linear dynamcs","title":"Faster convergence for linear dynamcs","text":"The optimized version simply has negative diagonal terms.","category":"page"},{"location":"02_dynamics/","page":"Faster convergence for linear dynamcs","title":"Faster convergence for linear dynamcs","text":"This may appear a bit trivial. In the next part, I optimize a system excluding the diagonal.","category":"page"},{"location":"02_dynamics/#Optimization-that-excludes-the-diagonal","page":"Faster convergence for linear dynamcs","title":"Optimization that excludes the diagonal","text":"","category":"section"},{"location":"02_dynamics/","page":"Faster convergence for linear dynamcs","title":"Faster convergence for linear dynamcs","text":"Here I consider a matrix that is stable, but produces a large nonlinear amplification. It is generated by the function:","category":"page"},{"location":"02_dynamics/","page":"Faster convergence for linear dynamcs","title":"Faster convergence for linear dynamcs","text":"function rand_nonnormal(n::Integer,(ud::Real)=1.01)\n  mat = randn(n,n) ./ sqrt(n)\n  @show SSA.spectral_abscissa(mat)\n  mat = mat -  (1.1*SSA.spectral_abscissa(mat))*I\n  sh = schur(mat)\n  upd = diagm(0=>fill(1.0,n),1=>fill(ud,n-1))\n  return sh.vectors*upd*sh.Schur*inv(upd)*sh.vectors'\nend","category":"page"},{"location":"02_dynamics/","page":"Faster convergence for linear dynamcs","title":"Faster convergence for linear dynamcs","text":"Let's make one and see how it looks like","category":"page"},{"location":"02_dynamics/","page":"Faster convergence for linear dynamcs","title":"Faster convergence for linear dynamcs","text":"A = rand_nonnormal(n,1.0)\nx0=randn(n)\ntimes,dyn_t,dyn_norms = run_linear_dyn(A,x0,30.,0.5)\nplot(times,dyn_norms;\n  leg=false,linewidth=3,color=:black,xlabel=\"time\",ylabel=\"norm(x(t))\")","category":"page"},{"location":"02_dynamics/","page":"Faster convergence for linear dynamcs","title":"Faster convergence for linear dynamcs","text":"The norm is initally amplified, and decreases slowly. This is due to the non-normality of matrix A.","category":"page"},{"location":"02_dynamics/#Objective-function-that-excludes-diagonal","page":"Faster convergence for linear dynamcs","title":"Objective function that excludes diagonal","text":"","category":"section"},{"location":"02_dynamics/","page":"Faster convergence for linear dynamcs","title":"Faster convergence for linear dynamcs","text":"function objective_and_grad_nodiag(x::Vector{R},grad::Union{Nothing,Vector{R}},\n    n::Integer,ssa_eps::R,alloc::SSA.SSAAlloc,A0::Matrix{R}) where R\n  mat=reshape(x,(n,n))\n  for i in 1:n\n    mat[i,i]=A0[i,i] # diagonal copied from original matrix\n  end\n  gradmat = isnothing(grad) ? nothing : similar(mat)\n  obj = SSA.ssa_simple!(mat,gradmat,alloc,ssa_eps)\n  if !isnothing(grad)\n    for i in 1:n\n      gradmat[i,i] = 0.0 # diagonal has zero gradient\n    end\n    for i in eachindex(gradmat)\n      grad[i]=gradmat[i] # copy the gradient\n    end\n  end\n  return obj\nend;\nnothing #hide","category":"page"},{"location":"02_dynamics/#Optimizer-and-optimization","page":"Faster convergence for linear dynamcs","title":"Optimizer and optimization","text":"","category":"section"},{"location":"02_dynamics/","page":"Faster convergence for linear dynamcs","title":"Faster convergence for linear dynamcs","text":"The only difference from before is using objective_and_grad_nodiag rather than objective_and_grad_simple","category":"page"},{"location":"02_dynamics/","page":"Faster convergence for linear dynamcs","title":"Faster convergence for linear dynamcs","text":"const ssa_eps=0.001\nconst alloc = SSA.SSAAlloc(n)\nconst y0 = A[:];\n\nfunction objfun!(F,G,y)\n  λ = 1.0/length(y) # regularizer weight\n  obj=objective_and_grad_nodiag(y,G,n,ssa_eps,alloc,A) # add the regularizer\n  ydiffs = y.-y0\n  obj += 0.5*λ*mapreduce(x->x^2,+,ydiffs)\n  if !isnothing(G)\n    @. G += λ*ydiffs # gradient of regularizer\n  end\n  return obj\nend\n\nopt_out = optimize(Optim.only_fg!(objfun!),A[:],BFGS(),Optim.Options(iterations=50))\ny_opt=Optim.minimizer(opt_out)\nA_opt = reshape(y_opt,(n,n));\nnothing #hide","category":"page"},{"location":"02_dynamics/","page":"Faster convergence for linear dynamcs","title":"Faster convergence for linear dynamcs","text":"Now the optimized matrix looks remarkably similar to ro the original one, as shown below.","category":"page"},{"location":"02_dynamics/","page":"Faster convergence for linear dynamcs","title":"Faster convergence for linear dynamcs","text":"heatmap(hcat(A,fill(NaN,n,10),A_opt);ratio=1,axis=nothing,ticks=nothing,border=:none)","category":"page"},{"location":"02_dynamics/","page":"Faster convergence for linear dynamcs","title":"Faster convergence for linear dynamcs","text":"Here I show the differences between A and A_textopt. (in case you don't believe they are different)","category":"page"},{"location":"02_dynamics/","page":"Faster convergence for linear dynamcs","title":"Faster convergence for linear dynamcs","text":"heatmap(A-A_opt;ratio=1,axis=nothing,ticks=nothing,border=:none)","category":"page"},{"location":"02_dynamics/","page":"Faster convergence for linear dynamcs","title":"Faster convergence for linear dynamcs","text":"Let's compare the time evolution of the norms","category":"page"},{"location":"02_dynamics/","page":"Faster convergence for linear dynamcs","title":"Faster convergence for linear dynamcs","text":"times,dyn_t_opt,dyn_norms_opt = run_linear_dyn(A_opt,x0,30.,0.5)\nplot(times,[dyn_norms dyn_norms_opt];\n      leg=:topright,linewidth=3,color=[:black :blue],\n      xlabel=\"time\",ylabel=\"norm(x(t))\", label=[\"before otpimization\" \"after optimization\"])","category":"page"},{"location":"02_dynamics/","page":"Faster convergence for linear dynamcs","title":"Faster convergence for linear dynamcs","text":"(Image: So much stability!)","category":"page"},{"location":"02_dynamics/#Extras","page":"Faster convergence for linear dynamcs","title":"Extras","text":"","category":"section"},{"location":"02_dynamics/","page":"Faster convergence for linear dynamcs","title":"Faster convergence for linear dynamcs","text":"In gradient based optimization with no automatic differentiation, it is always necessary to test the gradient of the objective function. The procedure is illustrated below.","category":"page"},{"location":"02_dynamics/","page":"Faster convergence for linear dynamcs","title":"Faster convergence for linear dynamcs","text":"using Calculus\nfunction test_gradient(myobjfun,y0)\n  grad_an = similar(y0)\n  _ = myobjfun(1.0,grad_an,y0) # compute gradient analytically\n  grad_num = Calculus.gradient(y->myobjfun(1.0,nothing,y),y0) # compute it numerically\n  return (grad_an,grad_num)\nend\n\n_ = let (x1,x2)=test_gradient(objfun!,randn(n^2))\n  scatter(x1,x2;ratio=1)\n  plot!(identity)\nend","category":"page"},{"location":"02_dynamics/","page":"Faster convergence for linear dynamcs","title":"Faster convergence for linear dynamcs","text":"","category":"page"},{"location":"02_dynamics/","page":"Faster convergence for linear dynamcs","title":"Faster convergence for linear dynamcs","text":"This page was generated using Literate.jl.","category":"page"},{"location":"01_show_ssa/","page":"Comparison between SSA and SA","title":"Comparison between SSA and SA","text":"EditURL = \"https://github.com/dylanfesta/SmoothedSpectralAbscissa.jl/blob/master/examples/01_show_ssa.jl\"","category":"page"},{"location":"01_show_ssa/#Comparison-between-SSA-and-SA","page":"Comparison between SSA and SA","title":"Comparison between SSA and SA","text":"","category":"section"},{"location":"01_show_ssa/","page":"Comparison between SSA and SA","title":"Comparison between SSA and SA","text":"In this example, I compare changes in the SA and in the SSA for a matrix that varies parametrically.","category":"page"},{"location":"01_show_ssa/#Initialization","page":"Comparison between SSA and SA","title":"Initialization","text":"","category":"section"},{"location":"01_show_ssa/","page":"Comparison between SSA and SA","title":"Comparison between SSA and SA","text":"using Plots,NamedColors\nusing LinearAlgebra\nusing Random\nusing SmoothedSpectralAbscissa ; const SSA=SmoothedSpectralAbscissa\nRandom.seed!(0);\nnothing #hide","category":"page"},{"location":"01_show_ssa/","page":"Comparison between SSA and SA","title":"Comparison between SSA and SA","text":"the function below generates a random non-normal matrix","category":"page"},{"location":"01_show_ssa/","page":"Comparison between SSA and SA","title":"Comparison between SSA and SA","text":"function rand_nonnormal(n::Integer,(ud::Real)=1.01)\n  mat = randn(n,n) ./ sqrt(n)\n  sh = schur(mat)\n  upd = diagm(0=>fill(1.0,n),1=>fill(ud,n-1))\n  return sh.vectors*upd*sh.Schur*inv(upd)*sh.vectors'\nend;\nnothing #hide","category":"page"},{"location":"01_show_ssa/#Example-matrix-generated-parametrically","page":"Comparison between SSA and SA","title":"Example matrix generated  parametrically","text":"","category":"section"},{"location":"01_show_ssa/","page":"Comparison between SSA and SA","title":"Comparison between SSA and SA","text":"n = 100\nmat1 = rand_nonnormal(n,0.8)\nmat2 = randn(n,n) ./ sqrt(n)\nmat(θ) = @. θ*mat1 + (1-θ)*mat2\n\nthetas = range(0.0,1.0;length=100);\nnothing #hide","category":"page"},{"location":"01_show_ssa/","page":"Comparison between SSA and SA","title":"Comparison between SSA and SA","text":"Now we look at the spectral abscissa as a function of the parameter","category":"page"},{"location":"01_show_ssa/#Set-ϵ-for-the-SSA","page":"Comparison between SSA and SA","title":"Set ϵ for the SSA","text":"","category":"section"},{"location":"01_show_ssa/","page":"Comparison between SSA and SA","title":"Comparison between SSA and SA","text":"ssa_eps_vals = [0.005,0.001,0.0005];\nnothing #hide","category":"page"},{"location":"01_show_ssa/#Compute-and-plot-the-SA","page":"Comparison between SSA and SA","title":"Compute and plot the SA","text":"","category":"section"},{"location":"01_show_ssa/","page":"Comparison between SSA and SA","title":"Comparison between SSA and SA","text":"sas = map(θ->SSA.spectral_abscissa(mat(θ)),thetas)\nplt=plot(thetas,sas ; color=:black, linewidth=3,lab=\"SA\",\n  xlabel=\"θ\",ylabel=\"SA value\",leg=:top)","category":"page"},{"location":"01_show_ssa/#Compute-and-plot-the-SSA","page":"Comparison between SSA and SA","title":"Compute and plot the SSA","text":"","category":"section"},{"location":"01_show_ssa/","page":"Comparison between SSA and SA","title":"Comparison between SSA and SA","text":"mycols=cgrad([colorant\"DarkGreen\",colorant\"orange\"])\nfor ϵ in ssa_eps_vals\n  ssas = map(θ->SSA.ssa_simple(mat(θ),ϵ),thetas)\n  plot!(plt,thetas,ssas ; linewidth=2.5 ,\n    lab=\"SSA $ϵ\",palette=mycols)\nend\nplot!(plt,xlabel=\"θ\",ylabel=\"SA/SSA value\",leg=:top)","category":"page"},{"location":"01_show_ssa/","page":"Comparison between SSA and SA","title":"Comparison between SSA and SA","text":"This plot shows that the SSA is a smoothed version of the SA, and coverges to it as as ϵ decreases. Moreover if we reduce the SSA parametrically, we find (approximately) a good minimum for the SA, as well.","category":"page"},{"location":"01_show_ssa/","page":"Comparison between SSA and SA","title":"Comparison between SSA and SA","text":"","category":"page"},{"location":"01_show_ssa/","page":"Comparison between SSA and SA","title":"Comparison between SSA and SA","text":"This page was generated using Literate.jl.","category":"page"},{"location":"#Smoothed-Spectral-Abscissa-(SSA)","page":"Home","title":"Smoothed Spectral Abscissa (SSA)","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"(Image: )","category":"page"},{"location":"","page":"Home","title":"Home","text":"This package computes the smoothed spectral abscissa (SSA) of square matrices, and the associated gradient, as described in:","category":"page"},{"location":"","page":"Home","title":"Home","text":"The Smoothed Spectral Abscissa for Robust Stability Optimization , J Vanbiervliet et al , 2009. DOI: 10.1137/070704034","category":"page"},{"location":"","page":"Home","title":"Home","text":"The SSA is a smooth upper bound to the spectral abscissa of a matrix, that is, the highest real part of the eigenvalues.","category":"page"},{"location":"","page":"Home","title":"Home","text":"The current version implements only the \"simplified\" version of SSA, i.e. the one where input and output transformations are identity matrices.","category":"page"},{"location":"#Definition-of-SSA","page":"Home","title":"Definition of SSA","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"Consider a square real matrix  A, that regulates a dynamical system as follows: textd mathbfxtextdt = Amathbfx. The stability of the dynamical system can be assessed using the following quantity:","category":"page"},{"location":"","page":"Home","title":"Home","text":"f(As) = int_0^infty textdt left expleft( left(A-sIright)t right)right^2","category":"page"},{"location":"","page":"Home","title":"Home","text":"where left M right^2 = texttraceleft(M M^top right).","category":"page"},{"location":"","page":"Home","title":"Home","text":"For a given varepsilon, the SSA can be denoted as tildealpha_varepsilon(A). And satisfies the following  equality:","category":"page"},{"location":"","page":"Home","title":"Home","text":"f(Atildealpha_varepsilon(A)) = frac1varepsilon","category":"page"},{"location":"","page":"Home","title":"Home","text":"It is an upper bound to the spectral abscissa of A. Therefore a reduction of SSA to values  0 guarantees dynamical stability.","category":"page"},{"location":"#Usage","page":"Home","title":"Usage","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"This module does not export functions in the global scope. It is therefore convenient to shorten the module name as follows:","category":"page"},{"location":"","page":"Home","title":"Home","text":"using SmoothedSpectralAbscissa ; const SSA=SmoothedSpectralAbscissa","category":"page"},{"location":"","page":"Home","title":"Home","text":"Therefore one can type SSA.foo in place of SmoothedSpectralAbscissa.foo.","category":"page"},{"location":"","page":"Home","title":"Home","text":"The functions below compute the SSA (and its gradient) for a matrix A.","category":"page"},{"location":"","page":"Home","title":"Home","text":"SSA.ssa_simple","category":"page"},{"location":"#SmoothedSpectralAbscissa.ssa_simple","page":"Home","title":"SmoothedSpectralAbscissa.ssa_simple","text":"ssa_simple(A,ssa_eps=nothing) -> ssa\n\nComputes the smoothed spectral abscissa (SSA) of matrix A (with identity input-output weighting matrices). If ssa_eps is not specified, the default value is set by default_eps_ssa(A)\n\nArguments\n\nA::Matrix{<:Real}: A square matrix\nssa_eps::Union{Nothing,R}=nothing The varepsilon parameter associated with the SSA computation. If not specified, it is set by default_eps_ssa(A)\n\nReturns\n\nssa<:Real: this is the tildealpha_varepsilon(A)\n\n\n\n\n\n","category":"function"},{"location":"","page":"Home","title":"Home","text":"SSA.ssa_simple_withgradient","category":"page"},{"location":"#SmoothedSpectralAbscissa.ssa_simple_withgradient","page":"Home","title":"SmoothedSpectralAbscissa.ssa_simple_withgradient","text":"  ssa_simple_withgradient(A,ssa_eps) -> (ssa,gradient_matrix)\n\nComputes the Smoothed Spectral Abscissa (SSA) of matrix A and also its gradient with respect to each element of A. If ssa_eps is not specified, the default value is set by default_eps_ssa(A)\n\nIn case you need to call the SSA iteratively (e.g. gradient-based optimization), please consider pre-allocating memory and using ssa_simple!(...).\n\n\n\n\n\n","category":"function"},{"location":"#Examples","page":"Home","title":"Examples","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"Comparison of SA and SSA\nSSA as objective to stabilize linear dynamics","category":"page"},{"location":"#Advanced-Interface","page":"Home","title":"Advanced Interface","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"When the SSA is used as optimization objective, it is convenient to use the advanced interface to avoid memory reallocation. The memory is pre-allocated in an object of type SSA.SSAlloc, which can then be used to call the SSA.ssa_simple(...) function multiple times.","category":"page"},{"location":"","page":"Home","title":"Home","text":"SSA.SSAAlloc","category":"page"},{"location":"#SmoothedSpectralAbscissa.SSAAlloc","page":"Home","title":"SmoothedSpectralAbscissa.SSAAlloc","text":"    SSAAlloc(n::Integer) -> SSAAlloc\n\nMemory allocation to avoid re-allocating space for repeated SSA computations.\n\nArguments\n\nn::Integer: the size of the matrix (or matrices) on which the SSA will be computed\n\n\n\n\n\n    SSAAlloc(A::Matrix) -> SSAAlloc\n\nMemory allocation to avoid re-allocating space for repeated SSA computations. This corresponds to SSAAlloc(size(A,1))\n\n\n\n\n\n","category":"type"},{"location":"","page":"Home","title":"Home","text":"Once the space is allocated, the SSA and its gradient can be computed by the following functions","category":"page"},{"location":"","page":"Home","title":"Home","text":"SSA.ssa_simple!","category":"page"},{"location":"#SmoothedSpectralAbscissa.ssa_simple!","page":"Home","title":"SmoothedSpectralAbscissa.ssa_simple!","text":"ssa_simple!(A,grad,PQ,ssa_eps=nothing) -> ssa\n\nComputes the smoothed spectral abscissa (SSA) of matrix A (with identity input-output weighting matrices) and its gradient (optionally). If ssa_eps is not specified, the default value is set by default_eps_ssa(A)\n\nArguments\n\nA::Matrix{<:Real}: A square matrix\ngrad::Union{Nothing,Matrix{<:Real}}: nothing  skips the gradient computation. To compute the gradient, provide a matrix grad of the same size as A. The matrix will be rewritten in-place with entries fracpartial tildealpha_varepsilon(A)partial A_ij\nPQ::SSAAlloc : Pre-allocated matrices for iterative computations. Can be generated by SSAAlloc(A)\nssa_eps::Union{Nothing,R}=nothing The varepsilon parameter associated with the SSA computation. If not specified, it is set by default_eps_ssa(A)\n\nReturns\n\nssa<:Real: this is tildealpha_varepsilon(A)\n\n\n\n\n\n","category":"function"},{"location":"","page":"Home","title":"Home","text":"SSA.ssa_simple_newton!","category":"page"},{"location":"#SmoothedSpectralAbscissa.ssa_simple_newton!","page":"Home","title":"SmoothedSpectralAbscissa.ssa_simple_newton!","text":"ssa_simple_newton!(A,grad,PQ,ssa_eps=nothing) -> ssa\n\nEquivalent to ssa_simple!( ... ), but uses the Netwon algorightm to compute the SSA. It might be slightly more efficient.\n\n\n\n\n\n","category":"function"},{"location":"#Index","page":"Home","title":"Index","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"","category":"page"}]
}
