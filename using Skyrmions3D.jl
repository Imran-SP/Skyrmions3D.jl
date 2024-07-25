using Skyrmions3D
using GLMakie
GLMakie.activate!()
Makie.inline!(false)

nuc = Skyrmion( [61, 61, 61], [0.2, 0.2, 0.2], mpi = 0.5, Fpi=100, ee=6.5, boundary_conditions="dirichlet")
overview(nuc)

p1(z) = z;
q1(z) = 1;

function f1(r, R)
    if r > R
        return 0
    else
        return π * (1 - r / R)
    end
end

#R = 1

#f1(r) = f1(r, R);

f2(r) = 4*atan(exp(-r));

make_rational_map!(nuc, p1, q1, f2)
Baryon(nuc)
Energy(nuc)

print_berger_grad_e2_star_at_31(nuc)

e2sgradient_flow!(nuc,tolerance=0.1,checks=100, dt = 0.000004)
#set_metric!(nuc,0.9)


#interactive_flow(nuc)

Energy(nuc)
gradient_flow!(nuc,tolerance=0.1,checks=100, dt = 0.00000004)

plot_baryon_density(nuc)