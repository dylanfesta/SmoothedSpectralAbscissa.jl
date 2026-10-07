@testset "Workspace constructors" begin
    A = [-1.0 0.3; -0.2 -2.0]
    for workspace in (SSA.Workspace(2), SSA.Workspace(A), SSA.Workspace(Float32.(A)))
        @test typeof(workspace) === SSA.Workspace{Matrix{Float64},Vector{ComplexF32}}
        buffers = [getfield(workspace, name) for name in fieldnames(typeof(workspace))]
        @test all(isconcretetype, fieldtypes(typeof(workspace)))
        @test all(buffer -> size(buffer) == (2, 2), buffers[1:12])
        @test size(workspace.A_eigvals) == (2,)
        @test all(i -> all(j -> buffers[i] !== buffers[j], 1:i-1), 2:13)
        @test SSA.PQ_init!(workspace, A) === nothing
        @test isapprox(SSA.spectral_abscissa(workspace), SSA.spectral_abscissa(A); atol=1e-6)
    end
    @test_throws ArgumentError SSA.Workspace(0)
    @test_throws ArgumentError SSA.Workspace(-1)
    @test_throws ArgumentError SSA.Workspace(zeros(0, 0))
    @test_throws DimensionMismatch SSA.Workspace(zeros(2, 3))
end

@testset "Unified SSA and workspace reuse" begin
    rng = MersenneTwister(42)
    matrices = ([-1.0 0.3; -0.2 -2.0], randn(rng, 2, 2) - 2I)
    workspace = SSA.Workspace(2)
    buffers = map(name -> getfield(workspace, name), fieldnames(typeof(workspace)))
    for method in (SSA.OptimOrder2, SSA.OptimNewton)
        for epsilon in (nothing, 0.2)
            for A in matrices
                original = copy(A)
                expected, expected_gradient = SSA.ssa_withgradient(A, epsilon)
                numerical_gradient = reshape(Calculus.gradient(
                    x -> SSA.ssa(copy(reshape(x, 2, 2)), epsilon), vec(copy(A))), 2, 2)
                @test isapprox(expected_gradient, numerical_gradient; atol=1e-6, rtol=1e-4)
                @test isapprox(SSA.ssa(A, epsilon; optim_method=method), expected; atol=1e-10)
                @test isapprox(SSA.ssa(A, epsilon; workspace=nothing,
                    optim_method=method), expected; atol=1e-10)
                @test isapprox(SSA.ssa(A, epsilon; workspace=workspace,
                    optim_method=method), expected; atol=1e-10)
                for storage in (nothing, workspace)
                    gradient = fill(NaN, size(A))
                    value = SSA.ssa(A, epsilon; workspace=storage, grad=gradient,
                        optim_method=method)
                    @test isapprox(value, expected; atol=1e-10)
                    @test isapprox(gradient, expected_gradient; atol=1e-10)
                    value, allocated_gradient = SSA.ssa_withgradient(A, epsilon;
                        workspace=storage, optim_method=method)
                    @test isapprox(value, expected; atol=1e-10)
                    @test isapprox(allocated_gradient, expected_gradient; atol=1e-10)
                    @test allocated_gradient !== A
                end
                @test A == original
                @test all(i -> getfield(workspace, fieldnames(typeof(workspace))[i]) === buffers[i], 1:13)
            end
        end
        value, gradient = SSA.ssa_withgradient([-2.0;;], 0.2; optim_method=method)
        @test isapprox(value, -1.9; atol=1e-12)
        @test isapprox(gradient, ones(1, 1); atol=1e-12)
    end
    @test !isdefined(SSA, :SSAAlloc)
    @test !isdefined(SSA, :ssa!)
    @test !isdefined(SSA, :ssa_simple!)
    @test_throws MethodError SSA.ssa(matrices[1], 0.2, I)
    @test_throws MethodError SSA.ssa_withgradient(matrices[1], 0.2, I)
end

@testset "SSA input validation" begin
    A = [-1.0 0.3; -0.2 -2.0]
    original = copy(A)
    workspace = SSA.Workspace(A)
    gradient = fill(-99.0, size(A))
    @test SSA.validate_ssa_inputs(A, gradient, workspace, 0.2) === nothing
    for compute in (SSA.ssa, SSA.ssa_withgradient)
        @test_throws DimensionMismatch compute(zeros(2, 3))
        @test_throws DimensionMismatch compute(zeros(2, 3); workspace=workspace)
        @test_throws ArgumentError compute(zeros(0, 0))
        @test_throws ArgumentError compute(zeros(0, 0); workspace=workspace)
        @test_throws DimensionMismatch compute(A; workspace=SSA.Workspace(3))
        for epsilon in (0.0, -0.1, Inf, -Inf, NaN)
            @test_throws ArgumentError compute(A, epsilon; workspace=workspace)
        end
    end
    @test_throws DimensionMismatch SSA.ssa(A; workspace=workspace, grad=zeros(3, 3))
    @test_throws ArgumentError SSA.ssa(A; workspace=workspace, grad=A)
    @test_throws ArgumentError SSA.ssa(A; workspace=workspace, grad=workspace.P)
    buffers = map(name -> getfield(workspace, name), fieldnames(typeof(workspace)))
    for i in 1:13
        bad_buffers = Base.setindex(buffers, zeros(3, 3), i)
        if i == 13
            bad_buffers = Base.setindex(buffers, zeros(ComplexF32, 3), i)
        end
        @test_throws DimensionMismatch SSA.ssa(A; workspace=SSA.Workspace(bad_buffers...), grad=gradient)
    end
    alias_buffers = Base.setindex(buffers, A, 1)
    @test_throws ArgumentError SSA.ssa(A; workspace=SSA.Workspace(alias_buffers...))
    @test A == original
    @test gradient == fill(-99.0, size(A))
end
