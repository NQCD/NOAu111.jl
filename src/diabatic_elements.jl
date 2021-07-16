
repulsive(A, α, r_cut) = JuLIP.@analytic r -> A * (exp(-α*r) - exp(-α*r_cut))
morse(F, γ, r₀) = JuLIP.@analytic r -> F*(1 - exp(-γ * (r - r₀)))^2
coupling(A2, A3, γ, r_cut) = JuLIP.@analytic r -> -A2*(1/(1+A3*exp(γ*r)) - 1/(1+A3*exp(γ*r_cut)))
const zero_potential = JuLIP.@analytic r -> 0

function repulsiveAuN(A, α, r₀, cosθ, r_cut)
    JuLIP.@analytic(
        r -> A*(exp(-2α*(r-r₀)) - exp(-2α*(r_cut-r₀)) - 2cosθ^2*(exp(-α*(r-r₀)) - exp(-α*(r_cut-r₀))))
    )
end
image(D, C, zimage) = JuLIP.@analytic z -> -D / sqrt(C^2 + (z-zimage)^2)

struct DiabaticElement{V,Z} <: JuLIP.SitePotential
    potentials::V
    zlist::Z
end

function DiabaticElement(potentials)
    zlist = JuLIP.Potentials.ZList([:Au, :N, :O])
    DiabaticElement(potentials, zlist)
end

JuLIP.cutoff(::DiabaticElement) = rcutoff

function H00()

    V00AuO = repulsive(A₀, α₀, rcutoff)
    V00AuN = repulsive(B₀, β₀, rcutoff)
    V00NO = morse(F₀, γ₀, r₀NO)

    potentials = DefaultDict{Tuple, JuLIP.AnalyticFunction}(zero_potential)
    potentials[1,3] = V00AuO
    potentials[1,2] = V00AuN
    potentials[2,3] = V00NO

    DiabaticElement(potentials)
end

function H01()

    V01AuO = coupling(A₂, A₃, γ₂, rcutoff)
    V01AuN = coupling(B₂, B₃, γ₃, rcutoff)

    potentials = DefaultDict{Tuple, JuLIP.AnalyticFunction}(zero_potential)
    potentials[1,3] = V01AuO
    potentials[1,2] = V01AuN

    DiabaticElement(potentials)
end

function H11()

    V11AuO = repulsive(A₁, α₁, rcutoff)
    V11NO = morse(F₁, γ₁, r₁NO)

    potentials = DefaultDict{Tuple, JuLIP.AnalyticFunction}(zero_potential)
    potentials[1,3] = V11AuO
    potentials[2,3] = V11NO

    DiabaticElement(potentials)
end

function JuLIP.evaluate!(tmp, V::DiabaticElement, Rs, Zs, z0)
    Es = 0.0
    i0 = JuLIP.z2i(V, z0)
    for (R, Z) in zip(Rs, Zs)
        i = JuLIP.z2i(V, Z)
        f = select_potential(V, i0, i)
        r = norm(R)
        Es += JuLIP.evaluate!(tmp, f, r)
    end
    return Es / 2
end

function JuLIP.evaluate_d!(dEs, tmp, V::DiabaticElement, Rs, Zs, z0)
    i0 = JuLIP.z2i(V, z0)
    for (j, (R, Z)) in enumerate(zip(Rs, Zs))
        i = JuLIP.z2i(V, Z)
        f = select_potential(V, i0, i)
        r = norm(R)
        R̂ = R / r
        dEs[j] = JuLIP.evaluate_d(f, r) * R̂ / 2
    end
    return dEs
end

function select_potential(V::DiabaticElement, i, j)

    sortkey(i, j) = i < j ? (i,j) : (j,i)

    key = sortkey(i, j)
    return V.potentials[key]
end
