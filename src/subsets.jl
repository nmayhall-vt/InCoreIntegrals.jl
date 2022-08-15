using InCoreIntegrals

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
        f = zeros(length(list),length(list))
        @tensor begin
            f[p,q] += viirs[p,q,r,s] * (da+db)[r,s]
            f[p,s] -= .5*viqri[p,q,r,s] * da[q,r]
            f[p,s] -= .5*viqri[p,q,r,s] * db[q,r]
        end
        ints_i.h1 .+= f
    end
    return ints_i
end


