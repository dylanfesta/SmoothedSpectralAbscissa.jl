using SmoothedSpectralAbscissa ; const SSA=SmoothedSpectralAbscissa
using LinearAlgebra, Random
using Roots,Calculus
using Test
Random.seed!(1)


"""
    lyap_simple!(A::AbstractMatrix)
returns the naive implementation of the Lyapunov objective function and its
derivative. Can be used for testing, plotting, etc
"""
function lyap_simple!(A::AbstractMatrix)
    function  myf(s::Float64)::Float64
        P = lyap(A - UniformScaling(s),one(A))
        tr(P)
    end
    function d_myf(s::Float64)::Float64
        sc,Id = UniformScaling(s), one(A)
        P =  lyap(A - sc,Id)
        Q = lyap(A' - sc,Id)
        -2.0tr(P*Q)
    end
    myf,d_myf
end

"""
    ssa_test(A::Matrix{Float64})

This is the easiest implementation
for testing and benchmarking.
"""
function ssa_test(A::Matrix{Float64} ; ssa_eps=nothing)
    myf,d_myf = lyap_simple!(A)
    ssa_eps = something(ssa_eps,SSA.default_eps_ssa(A))
    myg(s)= 1.0 / myf(s) - ssa_eps
    d_myg(s) = - d_myf(s) / myf(s)^2
    sa = SSA.spectral_abscissa(A)
    ssa_start = sa + 0.5*ssa_eps
    find_zero( (myg,d_myg),ssa_start, Roots.Newton())
end


function f_for_gradient(vmat)
    n = sqrt(length(vmat))
    @assert isinteger(n)
    n = Int64(n)
    mat = copy(reshape(vmat,(n,n)))
    return SSA.ssa(mat)
end

#

@testset "SA and SSA" begin
    # check SA
    n = 100
    mat = randn(n,n) + 1.4517I
    alloc = SSA.SSAAlloc(n)
    @test begin
        sa1 = maximum(real.(eigvals(mat)))
        SSA.PQ_init!(alloc,mat)
        sa2 = SSA.spectral_abscissa(alloc)
        sa3 = SSA.spectral_abscissa(mat)
        isapprox(sa1,sa2 ; rtol=1E-4) && isapprox(sa1,sa3 ; rtol=1E-4)
    end
    # now test SSA against simpler version
    ssa1 = ssa_test(mat)
    ssa2 = SSA.ssa(mat)
    @test isapprox(ssa1,ssa2 ; rtol=1E-4)
    # same , but with a different epsilon
    mat = randn(n,n) + 1.456I
    ε = 0.31313131
    ssa1 = ssa_test(mat; ssa_eps=ε)
    ssa2 = SSA.ssa(mat,ε)
    @test isapprox(ssa1,ssa2 ; rtol=1E-4)
end

@testset "SSA Vs Lyapunov eq" begin
    # when the epsilon of the ssa is the inverse of the energy integral
    # (computed by the lyap eq) , then the SSA is zero
    n = 130
    idm =diagm(0=>fill(1.0,n))
    matd = diagm(0=>fill(-0.05,n))
    fval=tr(lyap(matd,idm))
    ssa=SSA.ssa!(matd,nothing,SSA.SSAAlloc(n),inv(fval))
    @test isapprox(ssa,0.0;atol=1E-6)
    matrand=randn(n,n) ./ sqrt(n) - 1.1I
    fval=tr(lyap(matrand,idm))
    ssa = SSA.ssa!(matrand,nothing,SSA.SSAAlloc(n),inv(fval))
    @test isapprox(ssa,0.0;atol=1E-6)
    matrand2 = matrand + 1.234I
    ssa = SSA.ssa!(matrand2,nothing,SSA.SSAAlloc(n),inv(fval))
    @test isapprox(ssa,1.234;atol=1E-6)
end

@testset "SSA Newton method" begin
    n = 50
    matrand=randn(n,n) ./ sqrt(n)
    alloc=SSA.SSAAlloc(n)
    SSA.PQ_init!(alloc,matrand)
    sa_mat = SSA.spectral_abscissa(matrand)
    eps_ssa = SSA.default_eps_ssa(matrand)
    fungrad(s) = SSA.ssa_simple_obj_newton(s,alloc,eps_ssa,sa_mat)[1]
    stests=range(sa_mat+0.01,3sa_mat;length=100)
    grad_num = map(s->Calculus.gradient(fungrad,s),stests)
    grad_an = let  vv = map(s->SSA.ssa_simple_obj_newton(s,alloc,eps_ssa,sa_mat),stests)
        map(v-> v[1]/v[2],vv)
    end
    @test all(isapprox.(grad_num,grad_an;atol=1E-5))
    ssa = SSA.ssa_simple!(matrand,nothing,alloc)
    ssa_ntw =  SSA.ssa!(matrand,nothing,alloc;optim_method=SSA.OptimNewton)
    @test isapprox(ssa,ssa_ntw;atol=1E-5)
end


@testset "SSA gradient" begin
    n = 26
    mat = rand(n,n) + UniformScaling(2.222)
    @test begin
        ssa1 = ssa_test(mat)
        ssa2 = SSA.ssa(mat)
        isapprox(ssa1,ssa2 ; rtol=1E-4)
    end
    @test begin
        ssa,grad_an = SSA.ssa_withgradient(mat)
        grad_num = Calculus.gradient(f_for_gradient,mat[:])
        all(isapprox.(grad_an[:],grad_num ; rtol=1E-3) )
    end
end

@testset "SSA input-output weighting" begin
    A = [-1.0 0.0; 0.0 -2.0]
    A_original = copy(A)
    alloc = SSA.SSAAlloc(A)
    grad = similar(A)

    for ssa_eps in (nothing, 0.2)
        expected_value, expected_grad = SSA.ssa_withgradient(A, ssa_eps)
        @test isapprox(SSA.ssa(A, ssa_eps, I), expected_value; rtol=1E-12)
        value, gradient = SSA.ssa_withgradient(A, ssa_eps, I)
        @test isapprox(value, expected_value; rtol=1E-12)
        @test isapprox(gradient, expected_grad; rtol=1E-12)

        for method in (SSA.OptimOrder2, SSA.OptimNewton)
            value = SSA.ssa!(A, grad, alloc, ssa_eps;
                optim_method=method, input_output_weighting=I)
            @test isapprox(value, expected_value; rtol=1E-10)
            @test isapprox(grad, expected_grad; rtol=1E-10)
            value_without_gradient = SSA.ssa!(A, nothing, alloc, ssa_eps;
                optim_method=method, input_output_weighting=I)
            @test isapprox(value_without_gradient, expected_value; rtol=1E-10)
        end
    end

    unsupported_weights = (Matrix{Float64}(I, 2, 2), [2.0 0.0; 0.0 3.0],
        Diagonal([1.0, 1.0]), 2I, 0I, 1.0I)
    for weighting in unsupported_weights
        @test_throws r"Input-output weighting.*not implemented" SSA.ssa(A, 0.2, weighting)
        @test_throws r"Input-output weighting.*not implemented" SSA.ssa_withgradient(A, 0.2, weighting)
        for method in (SSA.OptimOrder2, SSA.OptimNewton)
            fill!(grad, -99.0)
            @test_throws r"Input-output weighting.*not implemented" SSA.ssa!(
                A, grad, alloc, 0.2; optim_method=method,
                input_output_weighting=weighting)
            @test grad == fill(-99.0, size(A))
        end
    end
    @test A == A_original
end
