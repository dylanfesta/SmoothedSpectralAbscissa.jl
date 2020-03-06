using Pkg
Pkg.activate(joinpath(@__DIR__,".."))
using LinearAlgebra
using BenchmarkTools


##

spectral_abscissa(A::AbstractMatrix) =  maximum(real.(eigvals(A)))

foo(v::SVector{4,Float64}) =  v .+ 1.0
foo(v::MVector{4,Float64}) =  v .+= 1.0
foo(v) = v .+= 1.0

mah = [0., 0. ,1. ,2]
boh = SVector{4}(mah)
buh = MVector{4}(mah)
##
# trace of matrix product
function trace_of_product(A::AbstractMatrix,B::AbstractMatrix)
    n=size(A,1)
    acc=0.0
    for i in 1:n, j in 1:n
        @inbounds @fastmath acc += A[i,j]*B[j,i]
    end
    return acc
end

function trace_of_product_meh(A::AbstractMatrix,B::AbstractMatrix)
    return tr(A*B)
end
# function trace_of_product(A::SMatrix{N,N,T},B::SMatrix{N,N,T}) where {T<:Real,N<:Integer}
#     return tr(A*B)
# end
function trace_of_product(A::SMatrix,B::SMatrix)
    return tr(A*B)
end


boh1 = rand(30,30)
boh2 = rand(30,30)
mah1 = SMatrix{30,30}(boh1)
mah2 = SMatrix{30,30}(boh2)

@btime trace_of_product($boh1,$boh2)
@btime trace_of_product_meh($boh1,$boh2)
@btime trace_of_product_meh($mah1,$mah2)
@btime trace_of_product($mah1,$mah2)

## Shur !

mah1 = SMatrix{100,100}(boh1)

@btime schur($boh1)
@time schur!(boh1)


function foo(M::Array{R,2}) where {N}
    @info N
end

gemm(A,B,C) = BLAS.gemm!('N','N',0.13,A,B,0.0,C)
gemmnot(A,B,C) = C .= 0.13 .*  (A*B)
gemmnot2(A,B,C) = C .= (A*B)
gemmnot3(A,B,C) = mul!(C,A,B,0.13, 0.0)


boh1 = rand(100,100)
boh2 = rand(100,100)
boh3 =  rand(100,100)

@btime gemmnot($boh1,$boh2,$boh3)
@btime gemmnot2($boh1,$boh2,$boh3)
@btime gemm($boh1,$boh2,$boh3)
@btime gemmnot3($boh1,$boh2,$boh3)

##

inv2(x) = 1 / x

@btime inv(1.333)
@btime inv2(1.333)

##

gem1(A,B,C) = BLAS.gemm!('T','N',-1.0,A,B,0.0,C)
gem2(A,B,C) = mul!(C,transpose(A),B,-1.0,0.0)


boh1 = rand(100,100)
boh2 = rand(100,100)
boh3 =  rand(100,100)

@btime gem1($boh1,$boh2,$boh3)
@btime gem2($boh1,$boh2,$boh3)

##

function subtract_diag1!(A,s)
    n=size(A,1)
    @simd for i in 1:n
        @inbounds A[i,i] -= s
    end
    return nothing
end
function subtract_diag2!(A,s)
    A[diagind(A,0)] .-=s
   return nothing
end
function subtract_diag3!(A,s)
    A -= s*I
end


boh1 = rand(100,100)
boh2 = rand(100,100)

@btime subtract_diag1!($boh1,0.001)
@btime subtract_diag2!($boh1,0.001)
@btime subtract_diag3!($boh1,0.001)

boh1
