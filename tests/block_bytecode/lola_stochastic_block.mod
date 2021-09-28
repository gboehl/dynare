/* Stochastic version of block decomposed LOLA model.
   Check that policy functions are the same as in non-block version. */

@#define deterministic = false
@#define block = true
@#define mfs = 0
@#include "lola_common.inc"

[~, state_reorder] = sort(oo_.dr.state_var);

ref = load(['lola_stochastic' filesep 'Output' filesep 'lola_stochastic_results.mat']);

[~, ref_state_reorder] = sort(ref.oo_.dr.state_var);

/* NB: With block, the rows of ghx and ghu are in declaration order (and not in
       DR-order as in non-block mode) */

if max(max(abs(oo_.dr.ghx(:, state_reorder) - ref.oo_.dr.ghx(ref.oo_.dr.inv_order_var, ref_state_reorder)))) > 3e-9
    error('Error in ghx')
end

if max(max(abs(oo_.dr.ghu - ref.oo_.dr.ghu(ref.oo_.dr.inv_order_var, :)))) > 4e-8
    error('Error in ghu')
end


