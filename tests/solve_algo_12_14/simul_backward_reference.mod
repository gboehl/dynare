// Reference simulation with “simul_backward” and default “solve_algo” value

@#include "simul_backward_common.inc"

s = simul_backward_model(initialconditions, @{simlength}, innovations);

save simul_backward_ref.mat s
