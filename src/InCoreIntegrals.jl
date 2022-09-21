module InCoreIntegrals

# using
using TensorOperations
using QCBase

include("type_InCoreInts.jl")
include("transformations.jl")
include("computations.jl")
include("subsets.jl")

# exports
export InCoreInts
#export orbital_rotation!
export orbital_rotation
export subset 
#export compute_energy
#export n_orb

end
