
# Smoothed Spectral Abscissa (SSA)

This package computes the smoothed spectral abscissa (SSA) of square matrices, and the associated gradient, as described in:

> Vanbiervliet, J. et al. (2009) “The Smoothed Spectral Abscissa for Robust Stability Optimization,” *SIAM Journal on Optimization*, 20(1), pp. 156–171. Available at: [https://doi.org/10.1137/070704034](https://doi.org/10.1137/070704034).

The SSA is a smooth upper bound to the spectral abscissa of a matrix, that is, the highest real part of the eigenvalues.

The current version implements only the "simplified" version of SSA, i.e. the one where input and output transformations are identity matrices.

## Definition of SSA

Consider a square real matrix  ``A``, that regulates a dynamical system as follows:
``\text{d} \mathbf{x}/\text{d}t = A\,\mathbf{x}``. The stability of the dynamical system
can be assessed using the following quantity:

```math
f(A,s) := \int_0^{\infty} \text{d}t \;\left\| \exp\left( \left(A-sI\right)t \right)\right\|^2
```

where ``\left\| M \right\|^2 := \text{trace}\left(M \,M^\top \right)``.

For a given ``\varepsilon``, the SSA can be denoted as ``\tilde{\alpha}_\varepsilon(A)``.
By definition, it satisfies the following  equality:

```math
f(A,\tilde{\alpha}_\varepsilon(A)) = \frac{1}{\varepsilon}
```

The SSA is an upper bound to the spectral abscissa (SA) of ``A``.
If matrix ``A`` is modified so that its SSA is below 0, the SA of ``A`` will also be
negative, which guarantees the stability of the associated linear dynamics.

## Usage

This module does not export functions in the global scope. It is therefore convenient to
shorten the module name as follows:

```julia
using SmoothedSpectralAbscissa ; const SSA=SmoothedSpectralAbscissa
```

It then becomes possible to use the shorter notation `SSA.foo` in place of `SmoothedSpectralAbscissa.foo`.

The functions below compute the SSA (and its gradient) for a matrix ``A``.

```@docs
SSA.ssa
```

```@docs
SSA.ssa_withgradient
```

## Examples

1. [**Comparison of SA and SSA**](generated/01_show_ssa.md)
2. [**Stability-optimized linear systems**](generated/02_dynamics.md)
3. [**Optimization of excitatory/inhibitory recurrent neural network**](generated/03_EI.md)

## Reusing working storage

Use the same `ssa` function for individual calls and repeated computations.
Omitting `workspace` allocates working storage internally; supplying it reuses
its matrices. For optimization, also preallocate the gradient:

```julia
A = [-1.0 0.3; -0.2 -2.0]
workspace = SSA.Workspace(A)
gradient = similar(A)
value = SSA.ssa(A, 0.2; workspace=workspace, grad=gradient)
```

`ssa` returns a scalar, preserves `A`, and overwrites the supplied workspace and
gradient despite having no `!` suffix. `ssa_withgradient` returns a tuple and
allocates a gradient matrix; it also accepts `workspace` and `optim_method`.
Workspaces can be reused with different matrices of the same size, but must not
be shared by concurrent computations. Input, gradient, and workspace buffers
must not alias each other.

```@docs
SSA.Workspace
```

The smoothing parameter defaults to `0.01 * 150 / size(A, 1)`.

```@docs
SSA.default_eps_ssa
SSA.PQ_init!
```

### Migrating the API

Replace `SSA.SSAAlloc(A)` with `SSA.Workspace(A)` and
`SSA.ssa!(A, gradient, alloc, epsilon)` with
`SSA.ssa(A, epsilon; workspace=workspace, grad=gradient)`.
The legacy `ssa_simple!` wrapper has been removed. Pass weighting as
`input_output_weighting=I`, rather than as a third positional argument.
Only identity weighting is currently supported.

## Index

```@index
```
