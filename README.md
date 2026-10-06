# Smoothed Spectral Abscissa

[![Documentation](https://img.shields.io/badge/docs-dev-blue.svg)](https://dylanfesta.github.io/SmoothedSpectralAbscissa.jl/dev/)
[![CI](https://github.com/dylanfesta/SmoothedSpectralAbscissa.jl/actions/workflows/CI.yml/badge.svg)](https://github.com/dylanfesta/SmoothedSpectralAbscissa.jl/actions/workflows/CI.yml)
[![License: CC0 1.0](https://img.shields.io/badge/license-CC0%201.0-blue.svg)](LICENSE)

This Julia package computes the smoothed spectral abscissa (SSA) of real square
matrices and its gradient. For repeated computations, `SSAAlloc` provides reusable
working matrices. The current implementation uses identity input/output weighting
matrices (the simplified formulation without projections).

Requires **Julia 1.10 or later**. The package has no plotting dependencies;
the executable documentation uses Makie and CairoMakie in a separate environment.

The algorithm is described in:

> The Smoothed Spectral Abscissa for Robust Stability Optimization, J. Vanbiervliet et al., 2009. [DOI: 10.1137/070704034](https://doi.org/10.1137/070704034)

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
