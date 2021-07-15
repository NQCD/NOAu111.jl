using NOAu111
using Test
using JuLIP
using PyCall
using ASE
using NonadiabaticDynamicsBase
using NonadiabaticModels
using LinearAlgebra: norm

build = pyimport("ase.build")

slab = build.fcc111("Au", size=(2,2,3), vacuum=10.0)
no = build.molecule("NO")
build.add_adsorbate(slab, no, 5, "ontop")

at = JuLIP.Atoms(ASEAtoms(slab))
rattle!(at, 0.5)

@testset "Components" begin
    @test JuLIP.Testing.fdtest(NOAu111.AuAu(), at, verbose=false)
    @test JuLIP.Testing.fdtest(NOAu111.H00(), at, verbose=false)
    @test JuLIP.Testing.fdtest(NOAu111.H01(), at, verbose=false)
    @test JuLIP.Testing.fdtest(NOAu111.H11(), at, verbose=false)
end

@testset "Final model" begin
    R = ang_to_au.(Matrix(hcat(at.X...)))
    c = PeriodicCell(ang_to_au.(cell(at)'))

    model = NOAu(chemical_symbols(at), c)

    @test JuLIP.Testing.fdtest(x->potential(model, x), x->derivative(model, x), R, verbose=false)
end
