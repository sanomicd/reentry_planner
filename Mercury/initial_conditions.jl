"""
Script derived from: https://youtu.be/7BA7iVTRyO4?si=SoczBFJ0xGUhcaz9
"""

include("../Science/constants.jl")
include("../Science/conversions.jl")


m = 1350 #kg
Cd = 1.5
Ve = 7500 #m/s
S = 2.8 #m^2
γe = -2.0*deg2rad()
ha = collect(range(0, he, 1000))

V = Ve.*exp.((1 ./(2 .*β)).*(ρs ./sin.(γe)).*(S.*Cd./m).*exp.(-β.*ha))

tout = collect(range(0, 10000, 10000000))
r0 = 100000
v0 = sqrt(2*G*m/r0)
xinit = [r0, v0]


tspan = tout #collect(range(0.0, 10000.0, ))
stateout = rk4(ballistic_derivatives, xinit, tspan)



