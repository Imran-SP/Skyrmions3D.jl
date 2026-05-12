using Skyrmions3D
using GLMakie
using CSV
using DataFrames
using Serialization
using LinearAlgebra
using Base.Threads
GLMakie.activate!()
Makie.inline!(false)


nuc = Skyrmion([120, 120, 120], [0.1, 0.1, 0.1], mpi = 0, Fpi = 130, ee = 6.5, boundary_conditions = "dirichlet")

p1(z) = z
q1(z) = 1
f1(r) = 4 * atan(exp(-r))

make_rational_map!(nuc, p1, q1, f1)

Berger_Isospin(nuc)

gradient_flow!(nuc,tolerance=0.1,checks=50,dt=0.0004)

Energy(nuc)
set_metric!(nuc4,1.1)

nuc4 = deserialize("l120_B4M0_metric1.0")
vMOI = compute_current(nuc4,label="vMOI")
evecs = eigvals(vMOI)


arrested_newton_flow!(nuc, tolerance = 0.01, checks = 50, dt=0.01)

plot_overview(nuc)


p4(z) = z^4 + 2.0*sqrt(3.0)*im*z^2 + 1.0;
q4(z) = z^4 - 2.0*sqrt(3.0)*im*z^2 + 1.0;
f4(r) = pi*exp( -(r.^3)./12.0 )

nuc4 =  Skyrmion([120, 120, 120], [0.1, 0.1, 0.1], mpi = 0, Fpi=130, ee=6.5, boundary_conditions="dirichlet")
overview(nuc4)

make_rational_map!(nuc4, p4, q4, f4)

arrested_newton_flow!(nuc4, tolerance = 0.01, checks = 10, dt=0.005)

nuc = deserialize("01ls_120lp_B4M0_t001_1.1")

results = DataFrame(Metric = Float64[], Isospin = Float64[], Energy = Float64[])

evalresults = DataFrame(Metric = Float64[], l1 = Float64[], l2 = Float64[], l3 = Float64[])


for metric in 1.1:0.1:7
    set_metric!(nuc4, metric)
    arrested_newton_flow!(nuc4, tolerance = 0.01, checks = 50, dt=0.01)
    energy = 12*pi*pi*Energy(nuc4)
    isospin = Berger_Isospin(nuc4) 
    itensor = compute_current(nuc4, label="vMOI")
    evals = eigvals(itensor)
    push!(evalresults, (Metric = metric, l1 = evals[1], l2 = evals[2], l3 = evals[3]))
    push!(results, (Metric = metric, Isospin = isospin, Energy = energy))
    filename = "01ls_120lp_B4M0_t001_$(metric)"
    serialize(filename, nuc4)
end

CSV.write("01ls_120lp_B4M0_t001.csv", results)
CSV.write("01ls_120lp_B4M0_t001_evals.csv", evalresults)


plot_overview(nuc)


