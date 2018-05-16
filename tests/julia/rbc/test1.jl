# Modification of the path (for packages). Should be done in ~/.juliarc.jl with a fixed path instead.
if isempty(findin([abspath("../../../julia")], LOAD_PATH))
    unshift!(LOAD_PATH, abspath("../../../julia"))
end

using Dynare
using SteadyState


# Compile the rbc.mod file -> produce a module with the model definition.
@dynare "rbc1.mod"

# First call to the steady state routine (analytical)
@time SteadyState.steady!(model_, oo_)

# First call to the steady state routine (analytical)
@time SteadyState.steady!(model_, oo_)

steadyState = deepcopy(oo_.steady_state)

y_init = copy(steadyState)
y_init[1] = 1.1*steadyState[1]
y_init[2] = 0.9*steadyState[2]

# First call to the steady state routine (numerical)
println("First call to the numerical steady state routine")
@time SteadyState.steady!(model_, oo_, y_init)

# Check results
@assert maximum(abs.(oo_.steady_state-steadyState))<1e-6

yinit = deepcopy(steadyState)
yinit[1] = 1.1*steadyState[1]
yinit[2] = 0.9*steadyState[2]

# Second call to the steady state routine (numerical)
println("Second call to the numerical steady state routine")
@time SteadyState.steady!(model_, oo_, yinit)

params = model_.params

# change alpha
println("Change α")
model_.params = deepcopy(params)
model_.params[4] = max(min(1.0, params[4]*1.1), 0.0)
@time ys = SteadyState.steady(model_, oo_) # Analytical steady state
@time SteadyState.steady!(model_, oo_, steadyState)
@assert maximum(abs.(oo_.steady_state-ys))<1e-6

# change delta
println("Change δ")
model_.params = deepcopy(params)
model_.params[6] = max(min(1.0, params[6]*1.1), 0.0)
@time ys = SteadyState.steady(model_, oo_) # Analytical steady state
@time SteadyState.steady!(model_, oo_, steadyState)
@assert maximum(abs.(oo_.steady_state-ys))<1e-6

# change beta
println("Change β")
model_.params = deepcopy(params)
model_.params[1] = max(min(1-1e-6, params[1]*0.99), 0.0)
@time ys = SteadyState.steady(model_, oo_) # Analytical steady state
@time SteadyState.steady!(model_, oo_, steadyState)
@assert maximum(abs.(oo_.steady_state-ys))<1e-6

# change tau
println("Change τ")
model_.params = deepcopy(params)
model_.params[3] = params[3]/1.1
@time ys = SteadyState.steady(model_, oo_)
@time SteadyState.steady!(model_, oo_, steadyState)
@assert maximum(abs.(oo_.steady_state-ys))<1e-6

# change Epsilon
println("Change ϵ")
model_.params = deepcopy(params)
model_.params[5] = params[5]*1.1
@time ys = SteadyState.steady(model_, oo_)
@time SteadyState.steady!(model_, oo_, steadyState)
@assert maximum(abs.(oo_.steady_state-ys))<1e-6
