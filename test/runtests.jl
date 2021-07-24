using NOAu111
using Test
using JuLIP
using PyCall
using ASE
using NonadiabaticDynamicsBase
using NonadiabaticModels
using Unitful, UnitfulAtomic

build = pyimport("ase.build")

slab = build.fcc111("Au", size=(2,4,4), vacuum=10.0)
no = build.molecule("NO")
build.add_adsorbate(slab, no, 5, "ontop")

at = JuLIP.Atoms(ASEAtoms(slab))
rattle!(at, 0.5)

@testset "H11" begin
    @test JuLIP.Testing.fdtest(NOAu111.AuNCoupling(), at, verbose=false)
    @test JuLIP.Testing.fdtest(NOAu111.H11DiabaticElement(), at, verbose=false)
    @test JuLIP.Testing.fdtest(NOAu111.ImagePotential(at), at, verbose=false)
end

@testset "Components" begin
    @test JuLIP.Testing.fdtest(NOAu111.AuAu(), at, verbose=false)
    @test JuLIP.Testing.fdtest(NOAu111.H00(), at, verbose=false)
    @test JuLIP.Testing.fdtest(NOAu111.H01(), at, verbose=false)
    @test JuLIP.Testing.fdtest(NOAu111.H11(at), at, verbose=false)
end

@testset "Final model" begin
    R = austrip.(Matrix(hcat(at.X...))u"Å")
    c = NonadiabaticDynamicsBase.PeriodicCell(austrip.(cell(at)' .* u"Å"))

    model = NOAu(chemical_symbols(at), c, R)

    @test JuLIP.Testing.fdtest(x->potential(model, x), x->derivative(model, x), R, verbose=false)
end
