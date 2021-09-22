module NOAu111

using JuLIP: JuLIP
using JuLIP.Potentials: z2i, ZList
using StaticArrays: SMatrix, SVector
using LinearAlgebra: norm, Hermitian, dot
using DataStructures: DefaultDict
using NonadiabaticDynamicsBase: au_to_ang, eV_per_ang_to_au, eV_to_au
using NonadiabaticModels: NonadiabaticModels, DiabaticModels

export NOAu

include("parameters.jl")
include("au_au.jl")
include("diabatic_elements.jl")
include("h11.jl")

struct NOAu{V1,V2,V3,V4,V5,A} <: DiabaticModels.DiabaticModel
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

NonadiabaticModels.nstates(::NOAu) = 2

NOAu(jatoms, Nindex, Oindex) = NOAu(2, H00(), H11(jatoms), H01(), AuAu(), image(D, C, zimage), jatoms, Nindex, Oindex)

function NOAu(symbols::AbstractVector{Symbol}, cell, R)
    jatoms = JuLIP.Atoms(
        X=zeros(3,length(symbols)),
        M=JuLIP.atomic_mass.(symbols),
        Z=JuLIP.atomic_number.(symbols),
        cell=au_to_ang.(cell.vectors'),
        pbc=cell.periodicity
        )

    JuLIP.set_positions!(jatoms, au_to_ang.(R))
    Nindex = findfirst(isequal(7), jatoms.Z)
    Oindex = findfirst(isequal(8), jatoms.Z)
    NOAu(jatoms, Nindex, Oindex)
end

function NonadiabaticModels.potential(model::NOAu, R::AbstractMatrix)

    JuLIP.set_positions!(model.atoms, au_to_ang.(R))
    Au = eV_to_au(JuLIP.energy(model.AuAu, model.atoms))
    V11 = eV_to_au(JuLIP.energy(model.H00, model.atoms)) + Au
    V22 = eV_to_au(JuLIP.energy(model.H11, model.atoms) + ϕ - Eₐ) + Au
    V12 = eV_to_au(JuLIP.energy(model.H01, model.atoms))

    return Hermitian(SMatrix{2,2}(V11, V12, V12, V22))
end

function NonadiabaticModels.derivative!(model::NOAu, D::AbstractMatrix{<:Hermitian}, R::AbstractMatrix)

    JuLIP.set_positions!(model.atoms, au_to_ang.(R))
    Au = -JuLIP.forces(model.AuAu, model.atoms)
    D11 = -JuLIP.forces(model.H00, model.atoms) + Au
    D22 = -JuLIP.forces(model.H11, model.atoms) + Au
    D12 = -JuLIP.forces(model.H01, model.atoms)

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

end
