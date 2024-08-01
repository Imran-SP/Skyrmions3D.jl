using Skyrmions3D
using BenchmarkTools
using GLMakie
GLMakie.activate!()
Makie.inline!(false)


nuc = Skyrmion([60, 60, 60], [0.2, 0.2, 0.2], mpi = 0, Fpi=100, ee=6.5, boundary_conditions="dirichlet")
overview(nuc)

p1(z) = z;
q1(z) = 1;

set_metric!(nuc,1)

f1(r) = 4*atan(exp(-r));

make_rational_map!(nuc, p1, q1, f1)
Baryon(nuc)
Energy(nuc)

gradient_flow!(nuc,tolerance=0.1, checks=50 ,dt=0.0004)

plot_baryon_density(nuc)

@benchmark Energy(nuc)

plot_overview(nuc)




p4(z) = z^4 + 2.0*sqrt(3.0)*im*z^2 + 1.0;
q4(z) = z^4 - 2.0*sqrt(3.0)*im*z^2 + 1.0;
f4(r) = pi*exp( -(r.^3)./12.0 )

nuc4 =  Skyrmion([60, 60, 60], [0.2, 0.2, 0.2], mpi = 0, Fpi=100, ee=6.5, boundary_conditions="dirichlet")
overview(nuc4)

make_rational_map!(nuc4, p4, q4, f4)
Baryon(nuc4)
Energy(nuc4)

gradient_flow!(nuc4,tolerance=0.1,checks=50,dt=0.0002)

set_metric!(nuc4,0.5)

interactive_flow(nuc4)

plot_overview(nuc4)