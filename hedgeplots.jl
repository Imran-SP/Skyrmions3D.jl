using Skyrmions3D
using GLMakie
GLMakie.activate!()
Makie.inline!(false)


nuc = Skyrmion([120, 120, 120], [0.1, 0.1, 0.1], mpi = 0.5, Fpi=100, ee=6.5, boundary_conditions="dirichlet")
overview(nuc)

p1(z) = z;
q1(z) = 1;


f1(r) = 4*atan(exp(-r));

make_rational_map!(nuc, p1, q1, f1)
Baryon(nuc)
Energy(nuc)

p,dp,ddp1,ddp2= getders_local_np(nuc,61,61,61)

get_berger_grad_e2_star(p,dp,ddp1)

nuc.pion_field[61,61,61,:]

print(dp)

plot_baryon_density(nuc)

#phi_0 = nuc.pion_field[:, :, :, 3]
#x_vals = 1:121
#phi_0_x = [phi_0[61, 61, x] for x in x_vals]


#fig_x = Figure()
#ax_x = Axis(fig_x[1,1], title="phi_3 vs z", xlabel="z", ylabel="phi_3")
#lines!(ax_x, x_vals, phi_0_x)
#fig_x