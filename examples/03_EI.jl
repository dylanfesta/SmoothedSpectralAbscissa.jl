# # Optimiziation of an excitatory/inhibitory (E/I) recurrent network

#=
In this example, I optimize an RNN network with E/I units.
The dynamics is expressed as follows:
```math
\tau\, \frac{\text{d} \mathbf{u}}{\text{d}t} = - \mathbf{u} + 
W f\left(\mathbf{u}\right) + \mathbf{h}
```
Where $W$ is a connection matrix with no autapses that preserves the separation
between E and I units (Dale's law). There is also a sparseness parameter for $W$,
to increase the difficulty (and as proof of concept).
Then there is a leak term $-\mathbf{u}$
and a constant external input $\mathbf{h}$. Finally the activation function
$f\left( \cdot \right)$ is a rectified-linear function.
=#

# ## Initialization
using Plots,NamedColors
using LinearAlgebra,Statistics
using Random
using SmoothedSpectralAbscissa ; const SSA=SmoothedSpectralAbscissa
using OrdinaryDiffEq # to run the dynamics
Random.seed!(0)

iofunction(x::Real) = max(0.0,x)

function run_rnn_dynamics(u0::Vector{R},W::Matrix{R},h::Vector{R},
    tmax::R,dt::R=0.01;verbose=false) where R
  f = function (du,u,p,t)
    mul!(du,W,iofunction.(u))
    @. du = du - u + h
    return du
  end
  prob = ODEProblem(f,u0,(0.,tmax))
  solv = solve(prob,Tsit5();verbose=verbose,saveat=dt)
  ret_u =hcat(solv.u...)
  ret_norms = mapslices(norm,ret_u;dims=1)[:]
  return solv.t,ret_u,ret_norms
end;

# ## Build the starting weight matrix

const sparseness = 0.5
const ne = 35
const ni = 15
const ntot = ne+ni
Wmask = hcat(fill(1.,ntot,ne),fill(-1.,ntot,ni)) # 1 for E , -1 for I, 0. for no connection
for i in 1:ntot,j in 1:ntot
  if i==j || rand()<sparseness
    Wmask[i,j]=0.
  end
end
heatmap(Wmask;
  ratio=1,seriescolor=:bwr,
  clims=(-1,1),cbar=nothing,axis=nothing,ticks=nothing,border=:none)
#=
The matrix `Wmask` shown above specifies the connectivy 
(presence of a conneciton, and wether the neuron is E or I)
The full weight matrix $W$ prior to optimization is defined as follows:
=#
W0 = 2.0 .* rand(ntot,ntot) .* Wmask;
#=
Even with the leaky term, the dynamics $\mathbf{u}(t)$ is very unstable.
=#
u0 = 3.0.*randn(ntot)
h = 0.1 .* rand(ntot)
times,dyn_t,dyn_norms = run_rnn_dynamics(u0,W0,h,10.0,0.05)
plot(times,dyn_norms;
  leg=false,linewidth=3,color=:black,xlabel="time",ylabel="norm(u(t))",
  yscale=:log10,label="before optimization")

# Note that the y-scale is exponential here.

# ## Objective function that excludes diagonal
using Optim
#=
To keep the signs of $W_{i,j}$ I make a reparametrization as follows.
```math
W_{i,j} = s_j \; \exp\left(\beta_{i,j}\right)
```
I then need to propagate the gradient of the SSA. I will also inclide a
2-norm regularizer on the weights
=#
function objective_and_grad_constraints(x::Vector{R},grad::Union{Nothing,Vector{R}},
    n::Integer,ssa_eps::R,alloc::SSA.SSAAlloc,W0::Matrix{R}) where R
  betas=reshape(x,(n,n))
  mat = @. exp(betas) * W0  # transform, apply constraints
  gradmat = isnothing(grad) ? nothing : similar(mat)
  obj = SSA.ssa!(mat,gradmat,alloc,ssa_eps)
  if !isnothing(grad)
    gradmat .*= mat # propagate gradient for constraint
    for i in eachindex(gradmat)
      grad[i]=gradmat[i]
    end
  end
  return obj
end;

# ## Optimizer and optimization

const ssa_eps=0.001
const alloc = SSA.SSAAlloc(ntot)
const y0 = map(w-> w!=0. ? log(abs(w)) : 0. ,W0[:])

function objfun!(F,G,y)
  λ = 3.0/length(y) # regularizer calibration
  obj=objective_and_grad_constraints(y,G,ntot,ssa_eps,alloc,W0) # add the regularizer
  obj += 0.5*λ*mapreduce(x->x^2,+,y)
  if !isnothing(G)
    @. G += λ*y # gradient of regularizer
  end
  return obj
end

opt_out = optimize(Optim.only_fg!(objfun!),y0,BFGS(),Optim.Options(iterations=1_000))
y_opt=Optim.minimizer(opt_out)
W_opt = Wmask .* exp.(reshape(y_opt,(ntot,ntot))); # convert from beta to weight

# Comparison between inital matrix and the optimized version. 
# Note the E/I separation and the spontaneous symmetry.
_ = let wboth = hcat(W0,fill(NaN,ntot,10),W_opt)
  _wex = extrema(W_opt)
  _zero_rel =abs(_wex[1])/(_wex[2]-_wex[1])
  _cgrad = cgrad([:blue,:white,:red],[0, _zero_rel ,1.])
  heatmap(wboth;ratio=1,seriescolor=_cgrad,clim=_wex,
    axis=nothing,ticks=nothing,border=:none)
end

# Now, I consider the norm of $\mathbf u(t)$ in the optimized system...
times,dyn_t_opt,dyn_norms_opt = run_rnn_dynamics(u0,W_opt,h,40.0,0.05)
plot(times,dyn_norms_opt;
      leg=:topright,linewidth=3,color=:blue,
      xlabel="time",ylabel="norm(x(t))", label="after optimization")

#=
STABLE !

There seems to be a large initial amplificaiton, but the activity settles to a stable point.
=#

# ## Extras

#=
In gradient based optimization with no automatic differentiation, it is
always necessary to test the gradient of the objective function.
The procedure is illustrated below.
=#

using Calculus
function test_gradient(myobjfun,y0)
  grad_an = similar(y0)
  _ = myobjfun(1.0,grad_an,y0) # compute gradient analytically
  grad_num = Calculus.gradient(y->myobjfun(1.0,nothing,y),y0) # compute it numerically
  for (k,y) in enumerate(y0)
    if y == 0.
      grad_num[k] = 0.
   end end
  return (grad_an,grad_num)
end

#=
To enable testing, set `do_the_test=true`
It is advised to reduce the size of the system.
=#

(x1,x2)=test_gradient(objfun!,randn(ntot^2))
scatter(x1,x2;ratio=1)
plot!(identity)

# Literate.markdown("examples/03_EI.jl","docs/src";documenter=true,repo_root_url="https://github.com/dylanfesta/SmoothedSpectralAbscissa.jl/blob/master") #src
