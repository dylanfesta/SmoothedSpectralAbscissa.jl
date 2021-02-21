# Smoothed Spectral Abscissa (SSA)
[![](https://img.shields.io/static/v1?logo=GitHub&label=.&message=SmoothedSpectralAbscissa.jl&color=blue)](https://github.com/dylanfesta/SmoothedSpectralAbscissa.jl)

This package computes the smoothed spectral abscissa (SSA) of square matrices, and the associated gradient, as described in:

> The Smoothed Spectral Abscissa for Robust Stability Optimization , J Vanbiervliet et al , 2009. [DOI: 10.1137/070704034](https://doi.org/10.1137/070704034)

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
And satisfies the following  equality:
```math
f(A,\tilde{\alpha}_\varepsilon(A)) = \frac{1}{\varepsilon}
```

It is an upper bound to the spectral abscissa of ``A``. Therefore a reduction of SSA to values
``< 0`` guarantees dynamical stability.

## Usage

This module does not export functions in the global scope. It is therefore convenient to
shorten the module name as follows:
```julia
using SmoothedSpectralAbscissa ; const SSA=SmoothedSpectralAbscissa
```
Therefore one can type `SSA.foo` in place of `SmoothedSpectralAbscissa.foo`.

The functions below compute the SSA (and its gradient) for a matrix ``A``.

```@docs
SSA.ssa_simple
```

```@docs
SSA.ssa_simple_withgradient
```


## Examples

1. [**Comparison of SA and SSA**](./01_show_ssa.md)
2. [**SSA as objective to stabilize linear dynamics**](./02_dynamics.md)


## Advanced Interface

When the SSA is used as optimization objective, it is convenient to use the advanced
interface to avoid memory reallocation. The memory is pre-allocated in an object of type
`SSA.SSAlloc`, which can then be used to call the `SSA.ssa_simple(...)` function multiple times.

```@docs
SSA.SSAAlloc
```

Once the space is allocated, the SSA and its gradient can be computed by the following
functions
```@docs
SSA.ssa_simple!
```
```@docs
SSA.ssa_simple_newton!
```


## Index

```@index
```
