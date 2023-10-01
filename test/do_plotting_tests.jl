
using Skyrmions3D
using Makie

a_skyrmion = Skyrmion(10,0.4)

@test_throws Exception plot_field(a_skyrmion)
@test_throws Exception plot_overview(a_skyrmion)
@test_throws Exception plot_baryon_density(a_skyrmion)

p(z) = z
q(z) = 1
f(r) = pi*exp( -(r.^3)./12.0 )

make_rational_map!(a_skyrmion, p, q)

plot_field(a_skyrmion)

@test typeof(plot_field(a_skyrmion)) == Makie.Figure
@test typeof(plot_overview(a_skyrmion)) == Makie.Figure
@test typeof(plot_baryon_density(a_skyrmion)) == Makie.Figure
@test typeof(plot_baryon_density(a_skyrmion, juggling=true)) == Makie.Figure

activate_CairoMakie()

@test typeof(plot_overview(a_skyrmion)) == Makie.Figure


