@#include "forward_looking_common.inc"

PSI=1.1;
TAU=2;
BETA=0.9;
KAPPA=0.6;

steady;
check;

estimated_params;
PSI,   1.1, NORMAL_PDF, 1.1, 0.01;
TAU,   2, NORMAL_PDF, 2, 0.01;
BETA,  0.9, NORMAL_PDF, 0.9, 0.01;
KAPPA, 0.6, NORMAL_PDF, 0.6, 0.01;
end;

% with this calibration and just a single observable the identification tests run into several edge case
% due to the fact that the order condition for dMoments fails and
% the minimal system cannot be computed
varobs x;
identification;
identification(advanced=1);
identification(prior_mc=10);
identification(prior_mc=10,advanced=1);
