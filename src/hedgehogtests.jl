using Skyrmions3D
using GLMakie
GLMakie.activate!()
Makie.inline!(false)


nuc = Skyrmion([60, 60, 60], [0.2, 0.2, 0.2], mpi = 0.5, Fpi=100, ee=6.5, boundary_conditions="dirichlet")
overview(nuc)

p1(z) = z;
q1(z) = 1;


f1(r) = 4*atan(exp(-r));

make_rational_map!(nuc, p1, q1, f1)
Baryon(nuc)
Energy(nuc)

p,dp,ddp1,ddp2= getders_local_np(nuc,30,30,30)

get_berger_grad_e2_star(p,dp,ddp1)

print(dp)


#set_metric!(nuc,0.93)
#e2sgradient_flow!(nuc,steps=10,tolerance=0.1,checks=1, dt = 0.000004)

print(nuc.pion_field[30,30,30,:])

interactive_flow(nuc)

Energy(nuc)
gradient_flow!(nuc,tolerance=0.1,checks=100, dt = 0.0004)

plot_baryon_density(nuc)



