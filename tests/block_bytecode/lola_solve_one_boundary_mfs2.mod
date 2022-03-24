// Tests option mfs=2 with block

@#define deterministic = true
@#define block = true
@#define mfs = 2
@#include "lola_common.inc"

mfs0=load(['lola_solve_one_boundary' filesep 'Output' filesep 'lola_solve_one_boundary_results']);

if max(max(oo_.endo_simul-mfs0.oo_.endo_simul)) > 20*options_.dynatol.x
   error('Inconsistency with mfs=0')
end
