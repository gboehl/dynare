% RBC model  with  irreversible  investment  constraint, implemented using MCP tag.

var k, y, L, c, A, a, mu, i;

varexo epsilon;

parameters beta, theta, tau, alpha, psi, delta, rho, Astar;

@#include "rbcii-calibration.inc"

model(use_dll);
  a = rho*a(-1) + epsilon;
  A = Astar*exp(a);
  y = A*(alpha*k(-1)^psi+(1-alpha)*L^psi)^(1/psi);
  k = y-c + (1-delta)*k(-1);
  (1-theta)/theta*c/(1-L) - (1-alpha)*(y/L)^(1-psi);
  (c^theta*(1-L)^(1-theta))^(1-tau)/c -mu =
  beta*(c(+1)^theta*(1-L(+1))^(1-theta))^(1-tau)/c(+1)
  *(alpha*(y(+1)/k)^(1-psi)+1-delta)+mu(+1)*(1-delta);
  i=y-c;
  mu = 0 ⟂ i > 0;
end;


steady_state_model;
  a = epsilon/(1-rho);
  A = Astar*exp(a);
  Output_per_unit_of_Capital=((1/beta-1+delta)/alpha)^(1/(1-psi));
  Consumption_per_unit_of_Capital=Output_per_unit_of_Capital-delta;
  Labour_per_unit_of_Capital=(((Output_per_unit_of_Capital/A)^psi-alpha)
  /(1-alpha))^(1/psi);
  Output_per_unit_of_Labour=Output_per_unit_of_Capital/Labour_per_unit_of_Capital;
  Consumption_per_unit_of_Labour=Consumption_per_unit_of_Capital
  /Labour_per_unit_of_Capital;
  % Compute steady state of the endogenous variables.
  L=1/(1+Consumption_per_unit_of_Labour/((1-alpha)*theta/(1-theta)
  *Output_per_unit_of_Labour^(1-psi)));
  c=Consumption_per_unit_of_Labour*L;
  k=L/Labour_per_unit_of_Capital;
  y=Output_per_unit_of_Capital*k;
  i=delta*k;
  mu=0;
end;

steady;

shocks;
  var epsilon;
  stderr 0.10;
end;

extended_path(periods=200,lmmcp);

if any(oo_.endo_simul(strmatch('i',M_.endo_names,'exact'),:)<-1e-6)
    error('lmmcp tag did not work.')
end

ds = dseries('rbcii-sim-data.mat');
if isoctave
    tolerance=5e-5;
else
    tolerance=1e-6;
end

if any(abs(transpose(oo_.endo_simul(strmatch('i',M_.endo_names,'exact'),:))-ds.Investment.data)>tolerance)
    error('Simulation with lmmcp returns different results.')
end

delete rbcii-sim-data.mat
