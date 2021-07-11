module NOAu111

import JuLIP
import JuLIP.Potentials: z2i, ZList
using StaticArrays
using LinearAlgebra

export H00
export H11
export H01

include("parameters.jl")
include("au_au.jl")

repulsive(A, α, r_cut) = JuLIP.@analytic r -> A * (exp(-α*r) - exp(-α*r_cut))
morse(F, γ, r₀) = JuLIP.@analytic r -> F*(1 - exp(-γ * (r - r₀)))^2
image(D, C, zimage) = JuLIP.@analytic z -> -D / sqrt(C^2 + (z-zimage)^2)
coupling(A2, A3, γ, r_cut) = JuLIP.@analytic r -> -A2*(1/(1+A3*exp(γ*r)) - 1/(1+A3*exp(γ*r_cut)))

abstract type ComponentSitePotential <: JuLIP.SitePotential end
abstract type CompoundPotential <: JuLIP.SitePotential end

CustomPotential = Union{ComponentSitePotential, CompoundPotential}

function JuLIP.energy(V::CompoundPotential, at::JuLIP.AbstractAtoms; kwargs...) 
   tmp = [JuLIP.alloc_temp(V, at) for i in 1:JuLIP.nthreads()]
   return JuLIP.energy!(tmp, V, at; kwargs...) + extra_potential(V, at)
end

JuLIP.cutoff(V::CustomPotential) = V.cutoff

sortkey(i, j) = i < j ? (i,j) : (j,i)

struct H00{V,Z,T} <: ComponentSitePotential
    potentials::V
    zlist::Z
    cutoff::T
end

function H00()
    cutoff = rcutoff

    V00AuO = repulsive(A₀, α₀, cutoff)
    V00AuN = repulsive(B₀, β₀, cutoff)
    V00NO = morse(F₀, γ₀, r₀NO)

    potentials = Dict((1,3)=>V00AuO, (1,2)=>V00AuN, (2,3)=>V00NO)
    zlist = JuLIP.Potentials.ZList([:Au, :N, :O])

    H00(potentials, zlist, cutoff)
end

function evaluate_potential(V::CustomPotential, i, j, R)
    key = sortkey(i, j)
    if haskey(V.potentials, key)
        return V.potentials[key](R)
    else
        return 0.0
    end
end

function JuLIP.evaluate!(tmp, V::CustomPotential, Rs, Zs, z0)
    Es = 0.0
    i0 = JuLIP.z2i(V, z0)
    for (R, Z) in zip(Rs, Zs)
        i = JuLIP.z2i(V, Z)
        Es += evaluate_potential(V, i0, i, R)
    end
    return Es / 2
end

struct H11{V,I,Z,T} <: CompoundPotential
    potentials::V
    V11image::I
    ϕ::T
    Eₐ::T
    zlist::Z
    cutoff::T
end

function H11()
    cutoff = rcutoff
    V11image = image(D, C, zimage)
    V11AuO = repulsive(A₁, α₁, cutoff)
    V11AuN = repulsive(B₁, β₁, cutoff)
    V11NO = morse(F₁, γ₁, r₁NO)

    potentials = Dict((1,3)=>V11AuO, (1,2)=>V11AuN, (2,3)=>V11NO)

    zlist = JuLIP.Potentials.ZList([:Au, :N, :O])
    H11(potentials, V11image, ϕ, Eₐ, zlist, cutoff)
end

function extra_potential(V::H11, at)
    com = (at[1]*at.M[1]+at[2]*at.M[2]) / (at.M[1]+at.M[2])
    zcom = com[3]
    return V.V11image(zcom) + V.ϕ - V.Eₐ
end

struct H01{V,Z,T} <: ComponentSitePotential
    potentials::V
    zlist::Z
    cutoff::T
end

function H01()
    cutoff = rcutoff

    V01AuO = coupling(A₂, A₃, γ₂, cutoff)
    V01AuN = coupling(B₂, B₃, γ₃, cutoff)

    potentials = Dict((1,3)=>V01AuO, (1,2)=>V01AuN)

    zlist = JuLIP.Potentials.ZList([:Au, :N, :O])
    H01(potentials, zlist, cutoff)
end

end
