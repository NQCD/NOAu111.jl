
export AuAu

struct AuAu <: JuLIP.SitePotential end

JuLIP.cutoff(::AuAu) = JuLIP.lattice_constant(:Au)

function JuLIP.evaluate!(tmp, ::AuAu, Rs, Zs, z0)
    Es = 0
    if z0 == 79 # Ensure current site is gold atom
        slice = (Zs .== 79) # Identify gold atoms
        AuRs = @view Rs[slice] # Slice gold atoms
        I = get_nearest_neighbours(AuRs)
        for R in @view AuRs[I]
            i = select_orientation_index(R)
            R -= direction_vectors[i]
            A = force_tensors[i]
            Es += R' * A * R
        end
    end
    return Es / 4
end

function JuLIP.evaluate_d!(dEs, tmp, ::AuAu, Rs, Zs, z0)
    if z0 == 79 # Ensure current site is gold atom
        slice = (Zs .== 79) # Identify gold atoms
        AuRs = @view Rs[slice] # Slice gold atoms
        I = get_nearest_neighbours(AuRs)
        for i in 1:length(dEs)
            dEs[i] = zero(dEs[i])
        end
        for (j, R) in zip(I, @view AuRs[I])
            i = select_orientation_index(R)
            R -= direction_vectors[i]
            A = force_tensors[i]
            dEs[j] = A*R / 2
        end
    end
    return dEs
end

function get_nearest_neighbours(AuRs)
    if length(AuRs) >= 12
        I = partialsortperm(AuRs, 1:12, by=x->norm(x))
    else
        I = partialsortperm(AuRs, 1:9, by=x->norm(x))
    end
    return I
end

"""
Identifies which force tensor is to be used based on symmetry considerations.

We calculate the dot product between the current vector and all possible basis vectors.
The largest dot product identifies the correct tensor.
"""
function select_orientation_index(vec)
    _, index = findmax(dot.(Ref(vec), direction_vectors))
    return index
end
