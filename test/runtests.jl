using NOAu111
using Test
using JuLIP
using PyCall
using ASE
using NonadiabaticDynamicsBase
using NonadiabaticModels
using LinearAlgebra
using FiniteDiff

ase = pyimport("ase")

slab = ase.build.fcc111("Au", size=(2,2,3), vacuum=10.0)
no = ase.build.molecule("NO")
ase.build.add_adsorbate(slab, no, 5, "ontop")

at = JuLIP.Atoms(ASEAtoms(slab))
rattle!(at, 0.5)

@testset "Components" begin
    @test JuLIP.Testing.fdtest(AuAu(), at, verbose=false)
    @test JuLIP.Testing.fdtest(H00(), at, verbose=false)
    @test JuLIP.Testing.fdtest(H01(), at, verbose=false)
    @test JuLIP.Testing.fdtest(H11(), at, verbose=false)
end

@testset "Final model" begin
    R = ang_to_au.(Matrix(hcat(at.X...)))
    atoms = NonadiabaticDynamicsBase.Atoms(chemical_symbols(at))
    c = PeriodicCell(ang_to_au.(cell(at)'))

    model = NOAu(atoms, c)

    @test JuLIP.Testing.fdtest(x->potential(model, x), x->derivative(model, x), R, verbose=false)
end
