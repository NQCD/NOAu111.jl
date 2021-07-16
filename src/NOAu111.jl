module NOAu111

import JuLIP
import JuLIP.Potentials: z2i, ZList
using StaticArrays
using LinearAlgebra: norm, Hermitian, dot
using DataStructures: DefaultDict
using NonadiabaticDynamicsBase
using NonadiabaticModels

export NOAu

include("parameters.jl")
include("au_au.jl")
include("diabatic_elements.jl")

struct NOAu{V1,V2,V3,V4,V5,A} <: DiabaticModel
    n_states::Int
    H00::V1
    H11::V2
    H01::V3
    AuAu::V4
    image_potential::V5
    atoms::A
    Nindex::Int
    Oindex::Int
end

NOAu(jatoms, Nindex, Oindex) = NOAu(2, H00(), H11(), H01(), AuAu(), image(D, C, zimage), jatoms, Nindex, Oindex)

function NOAu(symbols::AbstractVector{Symbol}, cell)
    jatoms = JuLIP.Atoms(
        X=zeros(3,length(symbols)),
        M=JuLIP.atomic_mass.(symbols),
        Z=JuLIP.atomic_number.(symbols),
        cell=au_to_ang.(cell.vectors'),
        pbc=cell.periodicity
        )
    Nindex = findfirst(isequal(7), jatoms.Z)
    Oindex = findfirst(isequal(8), jatoms.Z)
    NOAu(jatoms, Nindex, Oindex)
end

function NonadiabaticModels.potential(model::NOAu, R::AbstractMatrix)

    JuLIP.set_positions!(model.atoms, au_to_ang.(R))
    update_angular_potential!(model, model.atoms)
    Au = eV_to_au(JuLIP.energy(model.AuAu, model.atoms))
    V11 = eV_to_au(JuLIP.energy(model.H00, model.atoms)) + Au
    V22 = eV_to_au(JuLIP.energy(model.H11, model.atoms) + ϕ - Eₐ) + Au
    V12 = eV_to_au(JuLIP.energy(model.H01, model.atoms))

    V22 += eV_to_au(evaluate_image_potential(model))

    return Hermitian(SMatrix{2,2}(V11, V12, V12, V22))
end

function NonadiabaticModels.derivative!(model::NOAu, D::AbstractMatrix{<:Hermitian}, R::AbstractMatrix)

    JuLIP.set_positions!(model.atoms, au_to_ang.(R))
    update_angular_potential!(model, model.atoms)
    Au = -JuLIP.forces(model.AuAu, model.atoms)
    D11 = -JuLIP.forces(model.H00, model.atoms) + Au
    D22 = -JuLIP.forces(model.H11, model.atoms) + Au
    D12 = -JuLIP.forces(model.H01, model.atoms)

    Oderiv, Nderiv = evaluate_image_derivative(model)
    D22[model.Nindex] += SVector{3}(0, 0, Nderiv)
    D22[model.Oindex] += SVector{3}(0, 0, Oderiv)

    for i=1:length(model.atoms)
        for j=1:3
            d11 = eV_per_ang_to_au(D11[i][j])
            d12 = eV_per_ang_to_au(D12[i][j])
            d22 = eV_per_ang_to_au(D22[i][j])
            D[j,i] = Hermitian(SMatrix{2,2}(d11, d12, d12, d22))
        end
    end

    return D
end

function evaluate_image_potential(model::NOAu)
    at = model.atoms
    Nindex = model.Nindex
    Oindex = model.Oindex
    total_mass = at.M[Oindex] + at.M[Nindex]
    zcom = (at[Oindex][3]*at.M[Oindex]+at[Nindex][3]*at.M[Nindex]) / total_mass
    image = JuLIP.evaluate(model.image_potential, zcom)
    return image
end

function evaluate_image_derivative(model::NOAu)
    at = model.atoms
    Nindex = model.Nindex
    Oindex = model.Oindex
    total_mass = at.M[Oindex] + at.M[Nindex]
    zcom = (at[Oindex][3]*at.M[Oindex]+at[Nindex][3]*at.M[Nindex]) / total_mass
    image = JuLIP.evaluate_d(model.image_potential, zcom)
    Oderiv = image * at.M[Oindex] / (total_mass)
    Nderiv = image * at.M[Nindex] / (total_mass)
    return Oderiv, Nderiv
end

function update_angular_potential!(model::NOAu, at)
    O = at[model.Oindex]
    N = at[model.Nindex]
    cosθ = (O[3] - N[3]) / norm(O - N)
    model.H11.potentials[1,2] = repulsiveAuN(B₁, β₁, r₁AuN, cosθ, rcutoff)
end

end
