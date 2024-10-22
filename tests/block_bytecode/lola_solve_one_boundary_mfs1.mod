// Tests option mfs=1 with block

@#define block = true
@#define mfs = 1
@#include "lola_common.inc"

mfs0=load(['lola_solve_one_boundary' filesep 'Output' filesep 'lola_solve_one_boundary_results']);

if max(max(oo_.endo_simul-mfs0.oo_.endo_simul)) > options_.dynatol.x
   error('Inconsistency with mfs=0')
end
