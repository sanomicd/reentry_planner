include("../Science/constants.jl")
#Grab mercury reentry constants
mercury_constants()

using MATLABDiffEq, ParameterizedFunctions

function ballistic_derivatives(t, state)
    
    #State Vec
    V = state[1]
    h = state[2]
    γ = state[3]
    
    #Density Model
    ρ = ρs*exp(-β*h)

    #Aerodynamic Model
    L=0
    D = 0.5*ρ*V^2*S*Cd

    #Gravity Model
    g = gs*(RE/(RE+h))^2
    r = RE+h

    #Dynamics
    dVdt = -D/m - g*sin(γ)
    dγdt = (L/m - (g-V^2/r)*cos(γ))/V
    dhdt = V*sin(γ)

    #Return statedot
    statedot = [dVdt, dhdt, dγdt]
end

function simp_derivative(t, state)
    statedot = 2*state
end