using InCoreIntegrals

"""
    -h0::Real                # constant energy shift
    -h1::Array{T,2}          # one electron integrals
    -h2::Array{T,4}          # two electron integrals (chemist's notation)

Type to hold a second quantized Hamiltonian coefficients in memory
"""
struct InCoreInts{T}
    h0::T
    h1::Array{T,2}
    h2::Array{T,4}
end

function InCoreInts(ints::InCoreInts{TT}, T::Type) where {TT}
    return InCoreInts{T}(T(ints.h0), Array{T,2}(ints.h1), Array{T,4}(ints.h2))
end




"""
    subset(ints::InCoreInts, list; rmd1a=nothing, rdm1b=nothing)

Extract a subset of integrals acting on orbitals in list, returned as `InCoreInts` type

# Arguments
- `ints::InCoreInts`: Integrals for full system 
- `list`: list of orbital indices in subset
- `rdm1a`: 1RDM for embedding α density to make CASCI hamiltonian
- `rdm1b`: 1RDM for embedding β density to make CASCI hamiltonian
"""
function subset(ints::InCoreInts{T}, list, rdm1a=nothing, rdm1b=nothing) where {T}
    ints_i = InCoreInts{T}(ints.h0, view(ints.h1,list,list), view(ints.h2,list,list,list,list))
    if rdm1b != nothing 
        if rdm1a == nothing
            throw(Exception)
        end
    end

    if rdm1a != nothing
        if rdm1b == nothing
            throw(Exception)
        end
        da = deepcopy(rdm1a)
        db = deepcopy(rdm1b)
        da[:,list] .= 0
        db[:,list] .= 0
        da[list,:] .= 0
        db[list,:] .= 0
        viirs = ints.h2[list, list,:,:]
        viqri = ints.h2[list, :, :, list]
        fa = zeros(length(list),length(list))
        fb = copy(fa)
        @tensor begin
            ints_i.h1[p,q] += viirs[p,q,r,s] * (da+db)[r,s]
            ints_i.h1[p,s] -= .5*viqri[p,q,r,s] * da[q,r]
            ints_i.h1[p,s] -= .5*viqri[p,q,r,s] * db[q,r]
        end
    end
    return ints_i
end





