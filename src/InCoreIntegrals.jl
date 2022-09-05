module InCoreIntegrals

# using
using TensorOperations

include("type_InCoreInts.jl")
include("transformations.jl")
include("computations.jl")
include("subsets.jl")

# exports
export InCoreInts
#export orbital_rotation!
export orbital_rotation
export compute_energy
export subset 
export n_orb

end
