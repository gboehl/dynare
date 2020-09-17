/* Verify that the “datafile” option of “perfect_foresight_setup” behaves as
   “initval_file” (see #1663) */

var c k;
varexo x;

parameters alph gam delt bet aa;
alph=0.5;
gam=0.5;
delt=0.02;
bet=0.05;
aa=0.5;


model;
c + k - aa*x*k(-1)^alph - (1-delt)*k(-1);
c^(-gam) - (1+bet)^(-1)*(aa*alph*x(+1)*k^(alph-1) + 1 - delt)*c(+1)^(-gam);
end;

initval;
x = 1;
k = ((delt+bet)/(1.0*aa*alph))^(1/(alph-1));
c = aa*k^alph-delt*k;
end;

steady;

perfect_foresight_setup(periods=200, datafile = ramst_initval_file_data_col_vec_mat);
if oo_.exo_simul(2) ~= 1.2
  error('datafile option problem with exogenous variable');
end
if oo_.endo_simul(2, 2) ~= 13
  error('datafile option problem with endogenous variable');
end

perfect_foresight_solver;
