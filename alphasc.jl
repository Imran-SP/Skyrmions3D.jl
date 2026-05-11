using Skyrmions3D
using GLMakie
using CSV
using DataFrames
using Serialization
using LinearAlgebra
using Base.Threads
GLMakie.activate!()
Makie.inline!(false)


nuc = Skyrmion([60, 60, 60], [0.2, 0.2, 0.2], mpi = 0.5, Fpi = 108, ee = 4.84, boundary_conditions = "dirichlet")

p1(z) = z
q1(z) = 1
f1(r) = 4 * atan(exp(-r))

make_rational_map!(nuc, p1, q1, f1)

gradient_flow!(nuc,tolerance=0.1,checks=50,dt=0.0004)

Energy(nuc)
set_metric!(nuc,1)

nuc4 = deserialize("l120_B4M0_metric1.0")
vMOI = compute_current(nuc,label="vMOI")
evecs = eigvecs(vMOI)


arrested_newton_flow!(nuc, tolerance = 0.01, checks = 50, dt=0.01)

plot_overview(nuc)


p4(z) = z^4 + 2.0*sqrt(3.0)*im*z^2 + 1.0;
q4(z) = z^4 - 2.0*sqrt(3.0)*im*z^2 + 1.0;
f4(r) = pi*exp( -(r.^3)./12.0 )

nuc4 =  Skyrmion([120, 120, 120], [0.1, 0.1, 0.1], mpi = 0, Fpi=130, ee=6.5, boundary_conditions="dirichlet")
overview(nuc4)

make_rational_map!(nuc4, p4, q4, f4)



nuc = deserialize("l120_B1M01.0")

results = DataFrame(Metric = Float64[], Isospin = Float64[], Energy = Float64[])

for metric in 1:0.1:10
    set_metric!(nuc4, metric)
    arrested_newton_flow!(nuc4, tolerance = 0.01, checks = 50, dt=0.005)
    energy = 12*pi*pi*Energy(nuc4)
    isospin = Berger_Isospin(nuc4) 
    push!(results, (Metric = metric, Isospin = isospin, Energy = energy))
    filename = "l120_B4M0  $(metric)"
    serialize(filename, nuc)
end

CSV.write("p2_110l_B1M0_t001.csv", results)


plot_overview(nuc)


