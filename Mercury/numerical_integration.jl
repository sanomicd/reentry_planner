function state_integration(derivatives, init_condition, tout)
    prob = ODEProblem(derivatives, init_condition, tout)
    sol = solve(prob)
    return sol.t, hcat(sol.u...)'
end


function rk4(f, y0, t)
    #f is the derivative function, y0 is the initial condition, and t is the t span
    n = length(t)
    y = zeros((n, length(y0)))
    y[1,:] = y0

    for i in 1:n-1
        h = t[i+1] - t[i]
        k1 = f(y[i, :], t[i])
        k2 = f(y[i,:] .+ k1 * h/2, t[i] + h/2)
        k3 = f(y[i,:] .+ k2 * h/2, t[i] + h/2)
        k4 = f(y[i,:] .+ k3 * h, t[i] + h)
        y[i+1,:] = y[i,:] .+ (h/6) .* (k1 + 2*k2 + 2*k3 + k4)
        #print("k1: "* string(k1)*", k2: "*string(k2)*", k3: "* string(k3) *", k4: "*string(k4))

        if isnan(y[i+1, 1])
            throw("Is NaN")
        end

        if (i+1)>n
            throw("Unexpected time expansion")
        end
    end

    return y
end