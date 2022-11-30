// Check the correctedness of perfect foresight simulation of a purely backward model with “solve_algo=14”

@#include "purely_backward_common.inc"

perfect_foresight_solver(solve_algo = 14);

ref = load('purely_backward_reference/Output/purely_backward_reference_results.mat');

if max(max(abs(oo_.endo_simul - ref.oo_.endo_simul))) > 5e-5
    error('Incorrect results for perfect foresight with solve_algo=14')
end


