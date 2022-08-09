# InCoreIntegrals

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://nmayhall-vt.github.io/InCoreIntegrals.jl/stable/)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://nmayhall-vt.github.io/InCoreIntegrals.jl/dev/)
[![Build Status](https://github.com/nmayhall-vt/InCoreIntegrals.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/nmayhall-vt/InCoreIntegrals.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![Coverage](https://codecov.io/gh/nmayhall-vt/InCoreIntegrals.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/nmayhall-vt/InCoreIntegrals.jl)

A simple package to store with the 1 and 2 electron integrals that appear in quantum chemistry calculations. 
Because the focus is on active-space type models where the number of molecular orbitals is not expected to grow beyond around 200, this just holds everything in memory. 

# Example
```julia
    ints_0 = npzread("data_h6/h0.npy");
    ints_1 = npzread("data_h6/h1.npy");
    ints_2 = npzread("data_h6/h2.npy");

    ints = InCoreInts(ints_0, ints_1, ints_2)

    # rotate orbitals by some random transformation
    A = rand(size(ints.h1)...)
    F = svd(A)
    U = F.U * F.Vt

    ints2 = orbital_rotation(ints, U);  
    orbital_rotation!(ints, U); # in place
```
In this example `ints_0` is a scalar, providing any energy constant (nuclear repulsion, core orbital energy, etc), 
`ints_1` is a matrix, $h_{pq}$, that contains the coefficients for a 1 body operator:

$$h_{pq}\hat{a}^\dagger_p\hat{a}_q.$$

and `ints_2` is a 4-index tensor, $g_{pqrs}$ that contains the coefficients for a 2-body operator:

$$g_{pqrs}\hat{a}^\dagger_p\hat{a}^\dagger_q\hat{a}_s\hat{a}_r.$$

The energy can be computed for a given 1 and 2 RDM with:
```julia
  compute_energy(ints::InCoreInts, rdm1, rdm2)
```
as

$$E = E0 + \sum_{pq}h_{pq}D_{pq} + .5 * \sum_{pqrs} g_{pqrs}\Gamma_{pqrs}$$
