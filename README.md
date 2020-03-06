# SmoothedSpectralAbscissa

<!-- [![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://dylanfesta.github.io/SmoothedSpectralAbscissa.jl/stable)
 -->
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://dylanfesta.github.io/SmoothedSpectralAbscissa.jl/dev)
[![Build Status](https://travis-ci.com/dylanfesta/SmoothedSpectralAbscissa.jl.svg?branch=master)](https://travis-ci.com/dylanfesta/SmoothedSpectralAbscissa.jl)
[![Codecov](https://codecov.io/gh/dylanfesta/SmoothedSpectralAbscissa.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/dylanfesta/SmoothedSpectralAbscissa.jl)


The Purpose of this package is computing the Smoothed Spectral Abscissa of square matrices, and the associated gradient. The computation can be optimized for recursion, so that each new calculation does not reallocate memory.

I implemented only the "simplified" version that does not have projections.

The algorithm is described in the following paper:

> The Smoothed Spectral Abscissa for Robust Stability Optimization , J Vanbiervliet et al , 2009. [DOI: 10.1137/070704034](https://doi.org/10.1137/070704034)

**WORK IN PROGRESS!**

[Documentation here](https://dylanfesta.github.io/SmoothedSpectralAbscissa.jl/dev)

## Citing

See `CITATION.bib` for the relevant reference(s).
