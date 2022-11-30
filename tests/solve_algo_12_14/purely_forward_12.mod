// Check the correctedness of perfect foresight simulation of a purely forward model with “solve_algo=12”

@#include "purely_forward_common.inc"

perfect_foresight_solver(solve_algo = 12);

ref = load('purely_forward_reference/Output/purely_forward_reference_results.mat');

if max(max(abs(oo_.endo_simul - ref.oo_.endo_simul))) > 2e-4
    error('Incorrect results for perfect foresight with solve_algo=12')
end


