# -*- coding: utf-8 -*-
using Pkg
Pkg.activate("GRIDAP", shared=true)


using Gridap

L = 3
domain = (0,L)
partition = (10,)
model = CartesianDiscreteModel(domain, partition);

model.grid.node_coords

labels = get_face_labeling(model)

add_tag_from_tags!(labels, "left", [1])
add_tag_from_tags!(labels, "right", [2])

labels

Ω = Triangulation(model)
Γ = Boundary(Ω, tags="right")

reffe = ReferenceFE(lagrangian, Float64, 1)

Vₕ = TestFESpace(Ω,reffe,dirichlet_tags="left")

Uₕ = TrialFESpace(Vₕ, 0.0) 

dΩ = Measure(Ω,2)
dΓ = Measure(Γ,2)

EA = 1.0e3
F = 10
q = -10;

a(u,w) = ∫(EA*(∇(w)⋅∇(u)))*dΩ
l(w) = ∫(q*w)dΩ + ∫(F*w)*dΓ;

op = AffineFEOperator(a,l,Uₕ,Vₕ);

uₕ = solve(op)

writevtk(Ω, "TEMP_solution", cellfields=["u"=>uₕ])

using GridapMakie, CairoMakie

set_theme!(theme_black())

plot(uₕ, figure = ( size = (500,200), ) )
