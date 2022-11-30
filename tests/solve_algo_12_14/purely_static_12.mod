// Check the correctedness of perfect foresight simulation of a purely static model with “solve_algo=12”

@#include "purely_static_common.inc"

perfect_foresight_solver(solve_algo = 12);

ref = load('purely_static_reference/Output/purely_static_reference_results.mat');

if max(max(abs(oo_.endo_simul - ref.oo_.endo_simul))) > 1e-16
    error('Incorrect results for perfect foresight with solve_algo=12')
end
