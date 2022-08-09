using InCoreIntegrals
using NPZ
using Test
using LinearAlgebra
using Random

@testset "InCoreInts" begin
    Random.seed!(2);
    
    ints_0 = npzread("data_h6/h0.npy");
    ints_1 = npzread("data_h6/h1.npy");
    ints_2 = npzread("data_h6/h2.npy");

    ints = InCoreInts(ints_0, ints_1, ints_2)

    @test isapprox(ints.h0, 7.045311130494659, atol=1e-16)

    # add tests for transformations
    A = rand(size(ints.h1)...)
    F = svd(A)
    U = F.U * F.Vt
    
    ints2 = orbital_rotation(ints, U);
    orbital_rotation!(ints, U);
    @test isapprox(ints.h1, ints2.h1)
    #display(tr(ints.h1))
    norb = size(ints.h1,1)
    #display(tr(reshape(ints.h2, (norb*norb, norb*norb))))
    
    @test isapprox(tr(ints.h1), -14.528372393862362)
    @test isapprox(tr(reshape(ints.h2, (norb*norb, norb*norb))), 5.572554617047113)

end
