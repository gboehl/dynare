/* Test for the initval_file() command. This file needs ramst_initval_file_data.m. It should give results similar to those of ramst.mod */

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

initval_file(filename = ramst_initval_file_data_row_vec_mat);
perfect_foresight_setup(periods=200);
if oo_.exo_simul(2) ~= 1.2
  error('initval_file problem with exogenous variable');
end
if oo_.endo_simul(2, 2) ~= 13
  error('initval_file option problem with endogenous variable');
end
perfect_foresight_solver;
if ~oo_.deterministic_simulation.status
   error('Perfect foresight simulation failed')
end

oo_.exo_simul = [];
oo_.endo_simul = [];

initval_file(filename = ramst_initval_file_data_col_vec_mat);

perfect_foresight_setup(periods=200);
if oo_.exo_simul(2) ~= 1.2
  error('initval_file problem with exogenous variable');
end
if oo_.endo_simul(2, 2) ~= 13
  error('initval_file problem with endogenous variable');
end
perfect_foresight_solver;
if ~oo_.deterministic_simulation.status
   error('Perfect foresight simulation failed')
end

if ispc()
    initval_file(filename = ramst_initval_file_excel);
    perfect_foresight_setup(periods=200);
    perfect_foresight_solver;
    if ~oo_.deterministic_simulation.status
       error('Perfect foresight simulation failed');
    end;
end

initval_file(datafile = 'ramst_data.m');
perfect_foresight_setup(periods = 200);
perfect_foresight_solver;
if ~oo_.deterministic_simulation.status
   error('Perfect foresight simulation failed')
end
oo_.exo_simul = [];
oo_.endo_simul = [];

initval_file(datafile = 'ramst_data.mat');
perfect_foresight_setup(periods = 200);
perfect_foresight_solver;
if ~oo_.deterministic_simulation.status
   error('Perfect foresight simulation failed')
end
oo_.exo_simul = [];
oo_.endo_simul = [];

initval_file(datafile = 'ramst_data.csv');
perfect_foresight_setup(periods = 200);
perfect_foresight_solver;
if ~oo_.deterministic_simulation.status
   error('Perfect foresight simulation failed')
end
oo_.exo_simul = [];
oo_.endo_simul = [];

initval_file(datafile = 'ramst_data.xlsx');
perfect_foresight_setup(periods = 200);
perfect_foresight_solver;
if ~oo_.deterministic_simulation.status
   error('Perfect foresight simulation failed')
end
oo_.exo_simul = [];
oo_.endo_simul = [];

if ispc;
  initval_file(datafile = 'ramst_data.cls');
  perfect_foresight_setup(periods = 200);
  perfect_foresight_solver;
  oo_.exo_simul = [];
  oo_.endo_simul = [];

  if ~oo_.deterministic_simulation.status
    error('Perfect foresight simulation failed')
  end;
end;

