// Scalar translation of deterministic_matrix.mod
// The declaration orders are the same as those of the expanded matrix version
var t X_1 X_2 X_3 Y_1 Y_2 Y_3 Z_1 Z_2 Z_3;
varexo E_1 E_2 E_3 U_1 U_2 U_3 V_1 V_2 V_3;
parameters P_1_1 P_1_2 P_1_3 P_2_1 P_2_2 P_2_3 P_3_1 P_3_2 P_3_3;

P_1_1 = .1;
P_1_2 = .2;
P_1_3 = -.3;
P_2_1 = -.1;
P_2_2 = .5;
P_2_3 = .1;
P_3_1 = .2;
P_3_2 = -.3;
P_3_3 = .9;

model;
  X_1 = P_1_1*X_1(-1) + P_1_2*X_2(-1) + P_1_3*X_3(-1) + E_1;
  X_2 = P_2_1*X_1(-1) + P_2_2*X_2(-1) + P_2_3*X_3(-1) + E_2;
  X_3 = P_3_1*X_1(-1) + P_3_2*X_2(-1) + P_3_3*X_3(-1) + E_3;
  Y_1 = P_1_1*Y_1(-1) + P_1_2*Y_2(-1) + P_1_3*Y_3(-1) + U_1;
  Y_2 = P_2_1*Y_1(-1) + P_2_2*Y_2(-1) + P_2_3*Y_3(-1) + U_2;
  Y_3 = P_3_1*Y_1(-1) + P_3_2*Y_2(-1) + P_3_3*Y_3(-1) + U_3;
  Z_1 = V_1;
  Z_2 = V_2;
  Z_3 = V_3;
  t = Y_1(-1) + Y_2(+1);
end;

shocks;
  var E_1;
  periods 1 2;
  values 0.5 0.2;

  var E_2;
  periods 1 2;
  values 0.6 0.3;

  var E_3;
  periods 1 2;
  values 0.7 0.4;

  var U_1;
  periods 3;
  values 0.8;

  var U_2;
  periods 3;
  values 0.7;

  var V_1;
  periods 1;
  values 0.1;

  var V_2;
  periods 1;
  values 0.1;

  var V_3;
  periods 1;
  values 0.1;
end;

perfect_foresight_setup(periods=50);
perfect_foresight_solver;

L = load('deterministic_matrix_results.mat');
if max(max(abs(L.oo_.endo_simul - oo_.endo_simul))) > 1e-12
  error('Failure in matrix expansion')
end
