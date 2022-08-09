using InCoreIntegrals

"""
    compute_energy(h0, h1, h2, rdm1, rdm2)

Given an energy shift `h0`, 1e integrals `h1`, and 2e ints `h2`
along with a 1rdm and 2rdm on the same space, return the energy
"""
function compute_energy(h0, h1, h2, rdm1, rdm2)
    e = h0
    e += sum(h1 .* rdm1)
    e += .5*sum(h2 .* rdm2)
    # @tensor begin
    # 	e  += .5 * (ints.h2[p,q,r,s] * rdm2[p,q,r,s])
    # end
    return e
end


"""
    compute_energy(ints::InCoreInts, rdm1, rdm2)

Return energy defined by `rdm1` and `rdm2`
"""
function compute_energy(ints::InCoreInts, rdm1, rdm2)
    e = ints.h0
    e += sum(ints.h1 .* rdm1)
    e += .5*sum(ints.h2 .* rdm2)
    # @tensor begin
    # 	e  += .5 * (ints.h2[p,q,r,s] * rdm2[p,q,r,s])
    # end
    return e
end
