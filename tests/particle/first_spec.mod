//File testing error message if initial state vector is not positive definite

@#include "first_spec_common.inc"

varobs q ca;

shocks;
var eeps = 0.04^2;
var nnu = 0.03^2;
var q = 0.01^2;
var ca = 0.01^2;
end;

stoch_simul(order=3,periods=200, irf=0);
send_endogenous_variables_to_workspace;
save('my_data.mat','q','ca');

estimation(datafile='my_data.mat',order=2,mode_compute=0,mh_replic=0,filter_algorithm=sis,nonlinear_filter_initialization=2
    ,cova_compute=0 %tell program that no covariance matrix was computed
);