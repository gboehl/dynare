// Check the correctedness of “simul_backward” with “solve_algo=12”

@#include "simul_backward_common.inc"

options_.solve_algo = 12;
s = simul_backward_model(initialconditions, @{simlength}, innovations);

ref = load('simul_backward_ref.mat');

if max(max(abs(s.data - ref.s.data))) > 5e-4
    error('Incorrect results for simul_backward with solve_algo=12')
end
