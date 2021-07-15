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
end

NOAu(jatoms) = NOAu(2, H00(), H11(), H01(), AuAu(), image(D, C, zimage), jatoms)

function NOAu(atoms::NonadiabaticDynamicsBase.Atoms{N,T}, cell) where {N,T}
    jatoms = JuLIP.Atoms{T}(
        X=zeros(3,length(atoms)),
        M=au_to_u.(atoms.masses),
        Z=JuLIP.AtomicNumber.(atoms.numbers),
        cell=au_to_ang.(cell.vectors'),
        pbc=cell.periodicity
        )
    NOAu(jatoms)
end

function NonadiabaticModels.potential(model::NOAu, R::AbstractMatrix)

    JuLIP.set_positions!(model.atoms, au_to_ang.(R))
    Au = eV_to_au(JuLIP.energy(model.AuAu, model.atoms))
    V11 = eV_to_au(JuLIP.energy(model.H00, model.atoms)) + Au
    V22 = eV_to_au(JuLIP.energy(model.H11, model.atoms) + ϕ - Eₐ) + Au
    V12 = eV_to_au(JuLIP.energy(model.H01, model.atoms))

    V22 += eV_to_au(evaluate_image_potential(model))

    return Hermitian(SMatrix{2,2}(V11, V12, V12, V22))
end

function NonadiabaticModels.derivative!(model::NOAu, D::AbstractMatrix{<:Hermitian}, R::AbstractMatrix)

    JuLIP.set_positions!(model.atoms, au_to_ang.(R))
    Au = -JuLIP.forces(model.AuAu, model.atoms)
    D11 = -JuLIP.forces(model.H00, model.atoms) + Au
    D22 = -JuLIP.forces(model.H11, model.atoms) + Au
    D12 = -JuLIP.forces(model.H01, model.atoms)

    Oderiv, Nderiv = evaluate_image_derivative(model)
    Nindex = findfirst(isequal(7), model.atoms.Z)
    Oindex = findfirst(isequal(8), model.atoms.Z)
    D22[Nindex] += SVector{3}(0, 0, Nderiv)
    D22[Oindex] += SVector{3}(0, 0, Oderiv)

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
    Nindex = findfirst(isequal(7), at.Z)
    Oindex = findfirst(isequal(8), at.Z)
    total_mass = at.M[Oindex] + at.M[Nindex]
    zcom = (at[Oindex][3]*at.M[Oindex]+at[Nindex][3]*at.M[Nindex]) / total_mass
    image = JuLIP.evaluate(model.image_potential, zcom)
    return image
end

function evaluate_image_derivative(model::NOAu)
    at = model.atoms
    Nindex = findfirst(isequal(7), at.Z)
    Oindex = findfirst(isequal(8), at.Z)
    total_mass = at.M[Oindex] + at.M[Nindex]
    zcom = (at[Oindex][3]*at.M[Oindex]+at[Nindex][3]*at.M[Nindex]) / total_mass
    image = JuLIP.evaluate_d(model.image_potential, zcom)
    Oderiv = image * at.M[Oindex] / (total_mass)
    Nderiv = image * at.M[Nindex] / (total_mass)
    return Oderiv, Nderiv
end

end
