# Smoothed Spectral Abscissa

[![Documentation](https://img.shields.io/badge/docs-dev-blue.svg)](https://dylanfesta.github.io/SmoothedSpectralAbscissa.jl/dev/)
[![CI](https://github.com/dylanfesta/SmoothedSpectralAbscissa.jl/actions/workflows/CI.yml/badge.svg)](https://github.com/dylanfesta/SmoothedSpectralAbscissa.jl/actions/workflows/CI.yml)
[![License: CC0 1.0](https://img.shields.io/badge/license-CC0%201.0-white.svg)](LICENSE)

This Julia package computes the smoothed spectral abscissa (SSA) of real square
matrices and its gradient. For repeated computations, `Workspace` provides reusable
working matrices. The current implementation uses identity input/output weighting
matrices (the simplified formulation without projections).

```julia
using SmoothedSpectralAbscissa
const SSA = SmoothedSpectralAbscissa
A = [-1.0 0.3; -0.2 -2.0]
value = SSA.ssa(A)
value, gradient = SSA.ssa_withgradient(A)

# Reuse storage across calls for matrices of the same size.
workspace = SSA.Workspace(A)
gradient = similar(A)
value = SSA.ssa(A, 0.2; workspace=workspace, grad=gradient)
```

`ssa` preserves `A` and overwrites supplied workspace and gradient storage.
Omitting the workspace allocates it internally. Computations currently require
`Matrix{Float64}` inputs.

This API replaces `SSAAlloc` with `Workspace` and `ssa!` with
`ssa(A, ssa_eps; workspace=workspace, grad=gradient)`. The legacy `ssa_simple!`
wrapper is removed. Input/output weighting is now a keyword argument on both
`ssa` and `ssa_withgradient`.

Requires **Julia 1.10 or later**. The package has no plotting dependencies;
the executable documentation uses Makie and CairoMakie in a separate environment.

The algorithm is described in:

> Vanbiervliet, J. et al. (2009) “The Smoothed Spectral Abscissa for Robust Stability Optimization,” *SIAM Journal on Optimization*, 20(1), pp. 156–171. Available at: [https://doi.org/10.1137/070704034](https://doi.org/10.1137/070704034).

**Work in progress.** [Documentation and usage](https://dylanfesta.github.io/SmoothedSpectralAbscissa.jl/dev/).

## Development

Run the package tests from the repository root:

```sh
julia --project=. -e 'using Pkg; Pkg.test()'
```

Build the documentation, including all executable examples:

```sh
julia --project=docs -e 'using Pkg; Pkg.develop(PackageSpec(path=".")); Pkg.instantiate()'
julia --project=docs docs/make.jl
```

Open `docs/build/index.html` to view the local build. Edit the Literate sources in
`examples/`; tutorial Markdown is regenerated during the build and is not committed.
Each source can also be run directly, for example:

```sh
julia --project=docs examples/01_show_ssa.jl
```
