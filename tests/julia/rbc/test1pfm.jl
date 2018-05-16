# Modification of the path (for packages). Should be done in ~/.juliarc.jl with a fixed path instead.
if isempty(findin([abspath("../../../julia")], LOAD_PATH))
    unshift!(LOAD_PATH, abspath("../../../julia"))
end

# Load packages
using Plots
using Dynare
using SteadyState
using PerfectForesightModelSolver

plotlyjs() # Choose a backend


# Compile the rbc.mod file -> produce a module with the model definition.

@dynare "rbc1.mod"


@time SteadyState.steady!(model_, oo_)

println(oo_.steady_state)

# Initialize paths for the endogenous variables
endogenousvariables = repmat(oo_.steady_state, 1, options_.pfmsolver.periods+2)
# Destroy part of the initial stock of physical capital.
endogenousvariables[1, 1] = .8*endogenousvariables[1, 1]

# Set path for the innovations (no shocks).
exogenousvariables = repmat(oo_.exo_steady_state', options_.pfmsolver.periods+2, 1)

# Simulate the transition to the steady state
@time simulate_perfect_foresight_model!(endogenousvariables, exogenousvariables, oo_.steady_state, model_, options_)

n = 200
dates = collect(0:n-1)
#plt[:figure](1)
plot(dates, vec(endogenousvariables[1,1:n]), color="black", linewidth=2.0, linestyle="-")

# Initialize paths for the endogenous variables
endogenousvariables = repmat(oo_.steady_state, 1, options_.pfmsolver.periods+2)
# Destroy part of the initial stock of physical capital...
endogenousvariables[1, 1] = .8*endogenousvariables[1, 1]
# ... and assume that TFP is initially above its steady state level.
endogenousvariables[6, 1] = .5

# Set path for the innovations (no shocks).
exogenousvariables = repmat(oo_.exo_steady_state', options_.pfmsolver.periods+2, 1)

# Simulate the transition to the steady state (we should have an hump shaped transition)
@time simulate_perfect_foresight_model!(endogenousvariables, exogenousvariables, oo_.steady_state, model_, options_)

n = 200
dates = collect(0:n-1)
plot(dates, vec(endogenousvariables[1,1:n]), color="red", linewidth=2.0, linestyle="-")

# Initialize paths for the endogenous variables
endogenousvariables = repmat(oo_.steady_state, 1, options_.pfmsolver.periods+2)

# Set path for the innovations (no shocks).
exogenousvariables = repmat(oo_.exo_steady_state', options_.pfmsolver.periods+2, 1)

# Assume positive expected TFP shock in period 10
exogenousvariables[10+1, 1] = 2

# Simulate the paths for the endogenous variables, given the expected shock
@time simulate_perfect_foresight_model!(endogenousvariables, exogenousvariables, oo_.steady_state, model_, options_)

# dates = collect(0:n-1)
# n = 200
# plt[:figure](2)
# subplot(221)
# title("Efficiency")
# plot(dates, vec(endogenousvariables[5,1:n]), color="black", linewidth=2.0, linestyle="-")
# subplot(223)
# title("Output")
# plot(dates, vec(endogenousvariables[2,1:n]), color="black", linewidth=2.0, linestyle="-")
# subplot(222)
# title("Consumption")
# plot(dates, vec(endogenousvariables[4,1:n]), color="black", linewidth=2.0, linestyle="-")
# subplot(224)
# title("Labour")
# plot(dates, vec(endogenousvariables[3,1:n]), color="black", linewidth=2.0, linestyle="-")
# suptitle("Expected positive expected shock")
