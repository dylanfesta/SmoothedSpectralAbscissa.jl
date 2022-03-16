# Smoothed Spectral Abscissa

| **Documentation**         | **Build Status**                                                      |
|:------------------------- |:--------------------------------------------------------------------- |
|  [![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://dylanfesta.github.io/SmoothedSpectralAbscissa.jl/dev) |  [![Build Status](https://github.com/dylanfesta/SmoothedSpectralAbscissa.jl/workflows/CI/badge.svg)](https://github.com/dylanfesta/SmoothedSpectralAbscissa.jl/actions)  [![Coverage](https://codecov.io/gh/dylanfesta/SmoothedSpectralAbscissa.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/dylanfesta/SmoothedSpectralAbscissa.jl) |

This package implements a memory-efficient algorithm to compute  the smoothed spectral abscissa (SSA) of square matrices, and the associated gradient. The computation can be optimized for recursion, so that each new calculation does not reallocate memory.

I implemented only the "simplified" version that does not have projections.

The algorithm is described in the following paper:

> The Smoothed Spectral Abscissa for Robust Stability Optimization , J Vanbiervliet et al , 2009. [DOI: 10.1137/070704034](https://doi.org/10.1137/070704034)

**WORK IN PROGRESS!**

[Documentation and usage](https://dylanfesta.github.io/SmoothedSpectralAbscissa.jl/dev)
