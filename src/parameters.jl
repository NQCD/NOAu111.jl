using Unitful

const energy_unit = u"kJ/mol" / Unitful.Na

# Energies converted to eV, distances in Angstrom, forces in eV/Angstrom.

const A₀ = ustrip(uconvert(u"eV", 457095 * energy_unit))
const α₀ = 3.7594
const B₀ = ustrip(uconvert(u"eV", 30707 * energy_unit))
const β₀ = 3.0082
const F₀ = ustrip(uconvert(u"eV", 638.5 * energy_unit))
const γ₀ = 2.743
const r₀NO = 1.15077

const A₁ = A₀
const α₁ = α₀
const B₁ = ustrip(uconvert(u"eV", 24.056 * energy_unit))
const β₁ = 1.9649
const r₁AuN = 2.3491
const D = ustrip(uconvert(u"eV", 347.22 * energy_unit))
const C = 1.2423
const zimage = 1.1536
const F₁ = ustrip(uconvert(u"eV", 495.98 * energy_unit))
const γ₁ = 2.4890
const r₁NO = 1.2904
const ϕ = ustrip(uconvert(u"eV", 511.37 * energy_unit))
const Eₐ = ustrip(uconvert(u"eV", 0.67540 * energy_unit))

const A₂ = ustrip(uconvert(u"eV", 11.842 * energy_unit))
const A₃ = ustrip(uconvert(u"eV", 0.0061803 * energy_unit))
const γ₂ = 1.3693
const B₂ = ustrip(uconvert(u"eV", 51.000 * energy_unit))
const B₃ = ustrip(uconvert(u"eV", 0.0047213 * energy_unit))
const γ₃ = 2.0194

const α = ustrip(uconvert(u"eV*Å^-2", -4.94u"N/m"))
const β = ustrip(uconvert(u"eV*Å^-2", 17.15u"N/m"))
const γ = ustrip(uconvert(u"eV*Å^-2", 19.40u"N/m"))

const rcutoff = 10.0

const D1 = SMatrix{3,3}([α 0 0;
                         0 β γ;
                         0 γ β])
const D2 = SMatrix{3,3}([α 0 0;
                         0 β -γ;
                         0 -γ β])
const D3 = SMatrix{3,3}([β 0 γ; 
                         0 α 0;
                         γ 0 β])
const D4 = SMatrix{3,3}([β 0 -γ;
                         0 α 0;
                         -γ 0 β])
const D5 = SMatrix{3,3}([β γ 0;
                         γ β 0;
                         0 0 α])
const D6 = SMatrix{3,3}([β -γ 0;
                         -γ β 0;
                         0 0 α])
const force_tensors = [D1, D2, D3, D4, D5, D6,
                       D1, D2, D3, D4, D5, D6]

const one = SVector{3}(0, 1, 1) * 2.04
const two = SVector{3}(0, 1, -1) * 2.04
const three = SVector{3}(1, 0, 1) * 2.04
const four = SVector{3}(-1, 0, 1) * 2.04
const five = SVector{3}(1, 1, 0) * 2.04
const six = SVector{3}(1, -1, 0) * 2.04
const direction_vectors = [one, two, three, four, five, six,
                          -one,-two,-three,-four,-five,-six]
