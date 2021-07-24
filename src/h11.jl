
using UnPack
using Zygote

struct H11{I,A,V} <: JuLIP.AbstractCalculator
    image::I
    AuN::A
    potential::V
end

function H11(jatoms)
    image = ImagePotential(jatoms)
    AuN = AuNCoupling()
    potential = H11DiabaticElement()
    H11(image, AuN, potential)
end

function H11DiabaticElement()
    V11AuO = repulsive(A₁, α₁, rcutoff)
    V11NO = morse(F₁, γ₁, r₁NO)

    potentials = DefaultDict{Tuple, JuLIP.AnalyticFunction}(zero_potential)
    potentials[1,3] = V11AuO
    potentials[2,3] = V11NO

    DiabaticElement(potentials)
end

function JuLIP.energy(calc::H11, at::JuLIP.AbstractAtoms; domain=1:length(at))
    e1 = JuLIP.energy(calc.image, at)
    e2 = JuLIP.energy(calc.AuN, at; domain=domain)
    e3 = JuLIP.energy(calc.potential, at; domain=domain)
    return e1 + e2 + e3
end

function JuLIP.forces(calc::H11, at::JuLIP.AbstractAtoms; kwargs...)
    f1 = JuLIP.forces(calc.image, at)
    f2 = JuLIP.forces(calc.AuN, at; kwargs...)
    f3 = JuLIP.forces(calc.potential, at; kwargs...)
    return f1 + f2 + f3
end

struct ImagePotential{V,T} <: JuLIP.AbstractCalculator
    V::V
    mO::T
    mN::T
    mtotal::T
    iN::Int
    iO::Int
    slab_plane::T
end

function ImagePotential(jatoms)
    iN = findfirst(isequal(7), jatoms.Z)
    iO = findfirst(isequal(8), jatoms.Z)
    mN = jatoms.M[iN]
    mO = jatoms.M[iO]
    mtotal = mO + mN

    gold_positions = jatoms.X[jatoms.Z .== 79]
    slab_plane = maximum(x->x[3], gold_positions)

    V = image(D, C, zimage)
    ImagePotential(V, mO, mN, mtotal, iN, iO, slab_plane)
end

JuLIP.energy(calc::ImagePotential, at::JuLIP.AbstractAtoms) = JuLIP.evaluate(calc, at.X)
JuLIP.forces(calc::ImagePotential, at::JuLIP.AbstractAtoms) = -JuLIP.evaluate_d(calc, at.X)

function JuLIP.evaluate(calc::ImagePotential, R)
    @unpack mO, mN, mtotal, iN, iO = calc
    zcom = (R[iO][3]*mO + R[iN][3]*mN) / mtotal - calc.slab_plane
    return JuLIP.evaluate(calc.V, zcom)
end

function JuLIP.evaluate_d(calc::ImagePotential, R)
    @unpack mO, mN, mtotal, iN, iO = calc
    zcom = (R[iO][3]*mO + R[iN][3]*mN) / mtotal - calc.slab_plane
    image = JuLIP.evaluate_d(calc.V, zcom)
    out = zero(R)
    out[iN] = SVector{3}(0, 0, image * mN / mtotal)
    out[iO] = SVector{3}(0, 0, image * mO / mtotal)
    return out
end

struct AuNCoupling{V,Z} <: JuLIP.SitePotential
    V::V
    zlist::Z
end

function AuNCoupling()
    V(cosθ, r) = B₁*(exp(-2β₁*(r-r₁AuN)) - exp(-2β₁*(rcutoff-r₁AuN))
                - 2cosθ^2*(exp(-β₁*(r-r₁AuN)) - exp(-β₁*(rcutoff-r₁AuN)))
            )
    zlist = JuLIP.Potentials.ZList([:Au, :N, :O])
    AuNCoupling(V, zlist)
end

function JuLIP.evaluate!(tmp, V::AuNCoupling, Rs, Zs, z0)
    Es = 0.0
    i0 = JuLIP.z2i(V, z0)
    if i0 == 2 # Evaluated only at N site
        cosθ = calculate_cosθ(Rs, Zs)
        for (R, Z) in zip(Rs, Zs)
            i = JuLIP.z2i(V, Z)
            if i == 1
                r = norm(R)
                Es += V.V(cosθ, r)
            end
        end
    end
    return Es
end

function calculate_cosθ(Rs, Zs)
    closest = find_closest(Rs, Zs, 8)
    NO_distance = Rs[closest]
    cosθ = NO_distance[3] / norm(NO_distance)
    return cosθ
end

"""
Find the nearest atom with `z == ztarget`.
This is used to make sure the NO molecule interacts only with itself,
not the periodic images.
"""
function find_closest(Rs, Zs, ztarget)
    target_indices = findall(isequal(ztarget), Zs)
    max_dist = 100
    closest = target_indices[1]
    for i in target_indices
        dist = norm(Rs[i])
        if dist < max_dist
            closest = i
            max_dist = dist
        end
    end
    return closest
end

"""
This is the slowest part of the code, a hand-coded derivative would be faster,
or possibly a different autodiff package.
"""
function JuLIP.evaluate_d!(dEs, tmp, V::AuNCoupling, Rs, Zs, z0)
    grad = Zygote.gradient(x -> JuLIP.evaluate!(tmp, V, x, Zs, z0), Rs)[1]
    if grad === nothing
        for i=1:length(dEs)
            dEs[i] = zero(dEs[i])
        end
    else
        for i=1:length(grad)
            if grad[i] === nothing
                dEs[i] = zero(dEs[i])
            else
                dEs[i] = grad[i]
            end
        end
    end
    return dEs
end

JuLIP.cutoff(::AuNCoupling) = rcutoff
