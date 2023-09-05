"""
Script derived from: https://youtu.be/7BA7iVTRyO4?si=SoczBFJ0xGUhcaz9
"""

include("../Science/constants.jl")
include("../Science/conversions.jl")
include("numerical_derivative.jl")
include("numerical_integration.jl")
mercury_constants()


import Pkg
if "Plots" in keys(Pkg.dependencies())
    print("Plots installed")
else
    Pkg.add("Plots")
end
using Plots

m = 1350 #kg
Cd = 1.5
Ve = 7500 #m/s
S = 2.8 #m^2
γe = -2.0*deg2rad()
ha = collect(range(0, he, 1000))

V = Ve.*exp.((1 ./(2 .*β)).*(ρs ./sin.(γe)).*(S.*Cd./m).*exp.(-β.*ha))

tout = collect(range(0, 10000, 100000))
r0 = 100000
v0 = sqrt(2*G*m/r0)
xinit = [Ve; he; γe]

# Pendulum test
#y0 = [pi-0.1; 0.0]
#t = collect(range(0, 10, 31))
#sol = rk4(pend_derivative, y0, t)

# Back to Mercury

sol = rk4(ballistic_derivatives, xinit, tout)

Vnum = sol[:,1]
hnum = sol[:,2]
γnum = sol[:,3]

plot(ha, V, label="Analytical", title="Alt vs Velocity")
plot!(hnum, Vnum, label="Numerical")

plot(tout, γnum.*1 ./deg2rad(), title="Flight Path Angle (deg)")