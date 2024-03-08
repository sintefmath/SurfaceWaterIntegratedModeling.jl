u0 = [1., 0.]
function fun2(du, u, p, t)
    du[2] = -u[1]
    du[1] = u[2]
end

tspan = (0.0, 10.0)
prob = ODEProblem(fun2, u0, tspan)

condition(u, t, integrator) = u[2] 
affect!(integrator) = terminate!(integrator)
cb = ContinuousCallback(condition, affect!)

sol = solve(prob, Tsit5(), callback=cb)


# ----------------------------------------------------------------------------


u0 = [0.1]
function dudt(du, u, p, t)
    du[1] = 2u[1]
end

tspan = (0.0, 4.0)
prob = ODEProblem(dudt, u0, tspan)

condition(u, t, integrator) = u[1] - 101
affect!(integrator) = terminate!(integrator)
cb = ContinuousCallback(condition, affect!)

#sol = solve(prob, Tsit5(), callback=cb)
sol = solve(prob, Tsit5(), callback=ContinuousCallback(condition, affect!))

