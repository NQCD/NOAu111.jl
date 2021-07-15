using NOAu111
using PyCall
using NonadiabaticDynamicsBase
using NonadiabaticModels
using NonadiabaticMolecularDynamics
using UnitfulAtomic, Unitful
using CairoMakie
using LinearAlgebra

build = pyimport("ase.build")

slab = build.fcc111("Au", size=(2,4,4), vacuum=10.0)
no = build.molecule("NO")
build.add_adsorbate(slab, no, 5, "hcp")
build.add_vacuum(slab, 30)

cell, atoms, R = extract_parameters_and_positions(slab)

model = NOAu(atoms.types, cell)

bond_length_ref = austrip(1.15077u"Å")
R[3,end-1] = 40
R[3,end] = 40 + bond_length_ref
reference = ustrip.(auconvert.(u"kJ", eigvals(potential(model, R))[1]) .* Unitful.Na)

function plot_figure_a!(ax, reference, bond, lims)
    bond_length = austrip(bond)
    zN = austrip.((collect(lims[1]:0.01:lims[2]) .+ 17.0668) .*u"Å" ) 
    traj = Matrix[]
    for z in zN
        positions = copy(R )
        positions[3,end-1] = z
        positions[3,end] = z + bond_length
        push!(traj, positions)
    end
    potentials = potential.(model, traj)
    potentials = [ustrip.(auconvert.(u"kJ", V) .* Unitful.Na) for V in potentials]
    ground_state = [eigvals(V)[1] for V in potentials] .- reference
    v11 = [V[1,1] for V in potentials] .- reference
    v12 = [V[1,2] for V in potentials]
    v22 = [V[2,2] for V in potentials] .- reference
    x = ustrip.(auconvert.(u"Å", zN)) .- 17.0668
    lines!(ax, x, v11, label="V11", color=:green, linestyle=:dash)
    lines!(ax, x, v22, label="V22", color=:red, linestyle=:dashdotdot)
    lines!(ax, x, -v12, label="V12", color=:purple, linestyle=:dashdot)
    lines!(ax, x, ground_state, label="E0", color=:blue)
    xlims!(ax, lims[1], lims[2])
    return ax
end

function plot_figure_e!(ax, reference, lims)
    z = austrip(1.6u"Å")
    rs = austrip.(collect(lims[1]:0.01:lims[2]).*u"Å")
    traj = Matrix[]
    positions = copy(R)
    positions[3, end-1] = z
    for r in rs
        p = copy(positions)
        p[3,end] = z + r
        push!(traj, p)
    end
    potentials = potential.(model, traj)
    potentials = [ustrip.(auconvert.(u"kJ", V) .* Unitful.Na) for V in potentials]
    ground_state = [eigvals(V)[1] for V in potentials] .- reference
    v11 = [V[1,1] for V in potentials] .- reference
    v12 = [V[1,2] for V in potentials]
    v22 = [V[2,2] for V in potentials] .- reference
    x = ustrip.(auconvert.(u"Å", rs))
    lines!(ax, x, v11, label="V11", color=:green, linestyle=:dash)
    lines!(ax, x, v22, label="V22", color=:red, linestyle=:dashdotdot)
    lines!(ax, x, -v12, label="V12", color=:purple, linestyle=:dashdot)
    lines!(ax, x, ground_state, label="E0", color=:blue)
    xlims!(ax, lims[1], lims[2])
    return ax
end

theme = Theme(linewidth=3)
set_theme!(theme)
f = Figure(resolution=(700, 1200), backgroundcolor = RGBf0(0.98, 0.98, 0.98))

ax1 = f[1,1] = Axis(f, yticks=LinearTicks(5))
label_a = f[1,1,TopRight()] = Label(f, "A", textsize=24)
label_a.padding = (-50, 0, -50, 0)

ax2 = f[2,1] = Axis(f)
label_b = f[2,1,TopRight()] = Label(f, "C", textsize=24)
label_b.padding = (-50, 0, -50, 0)

ax3 = f[3,1] = Axis(f)
label_c = f[3,1,TopRight()] = Label(f, "E", textsize=24)
label_c.padding = (-50, 0, -50, 0)

plot_figure_a!(ax1, reference, 1.191u"Å", (1.0, 2.5))
ylims!(ax1, -50, 400)

plot_figure_a!(ax2, reference, 1.6u"Å", (0.5, 3.0))
ylims!(ax2, -20, 800)

plot_figure_e!(ax3, reference, (0.8, 2.0))
ylims!(ax3, -100, 2000)

ax1.xlabel = "N atom height / Å"
ax2.ylabel = "Energy /kJ mol⁻¹"
ax2.xlabel = "N atom height / Å"
ax3.xlabel = "N-O bond length / Å"

CairoMakie.save("plot/fig2.png", f)
f
