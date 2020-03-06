# SmoothedSpectralAbscissa.jl

The Purpose of this package is computing the Smoothed Spectral Abscissa of square matrices, and the associated gradient. The computation can be optimized for recursion, so that each new calculation does not reallocate memory.

I implemented only the "simplified" version that does not have projections.

The algorithm is described in the following paper:

> The Smoothed Spectral Abscissa for Robust Stability Optimization , J Vanbiervliet et al , 2009. [DOI: 10.1137/070704034](https://doi.org/10.1137/070704034)



## TO-DOs

  + Add references
  + describe the algorithm   

--------------------

```@index
```


```@autodocs
Modules = [SmoothedSpectralAbscissa]
```
