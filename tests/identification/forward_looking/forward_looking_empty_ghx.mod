@#include "forward_looking_common.inc"

PSI=1.1;
TAU=2;
BETA=0.9;
KAPPA=0.6;

steady;
check;

estimated_params;
PSI,   1.1;
TAU,   2;
BETA,  0.9;
KAPPA, 0.6;
end;

varobs r x p;
identification; %this triggers sylvester3a with empty ghx

% as a side note, we have the true solution:
% [r;x;p] = TRUE_SOLUTION*[e_M;e_D;e_S] (ghx is empty)
A = [1 0 -PSI; 1/TAU 1 0; 0 -KAPPA 1];
TRUE_SOLUTION1 = inv(A);
TRUE_SOLUTION2 = 1/(KAPPA*PSI/TAU +1)*[1         KAPPA*PSI PSI;
                                      -1/TAU     1         -PSI/TAU;
                                      -KAPPA/TAU KAPPA     1];
% note that BETA drops out from the solution

if max(max(abs(TRUE_SOLUTION1 - oo_.dr.ghu))) > 1e-15
    error('Something wrong with perturbation');
end
if max(max(abs(TRUE_SOLUTION2 - oo_.dr.ghu))) > 1e-15
    error('Something wrong with perturbation');
end