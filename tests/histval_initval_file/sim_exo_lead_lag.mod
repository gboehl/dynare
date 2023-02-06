// Uses autonomous system from sim_base.mod, but adds separate system where exogenous variables have several leads and lags
// Lags and leads on exogenous variables are substituted out by auxiliary variables

var c cmav k z_backward z_forward;
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
  z_backward=0.1*1+0.9*z_backward(-1) + (x(-4) - 1);
  z_forward=0.2*1+0.8*z_forward(+1) + (x(+4) - 1);
  cmav = 0.2*(c(-2) + c(-1) + c + c(+1) + c(+2));
end;

initval;
  c = 1.2;
  cmav = 1.2;
  k = 12;
  x = 1; //set x(0), x(-1), x(-2), x(-3)
  z_backward = 1;
  z_forward = 1;
end;

shocks;
var x; //sets x(+2)
periods 2;
values 0.9;
end;

perfect_foresight_setup(periods=200);
perfect_foresight_solver;

if ~oo_.deterministic_simulation.status
   error('Perfect foresight simulation failed')
end


