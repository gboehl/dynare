// Scalar translation of stochastic_matrix.mod
// The declaration orders are the same as those of the expanded matrix version
var r X_1 X_2 X_3 Y_1 Y_2 Y_3 Z_1 Z_2 Z_3 T_1 T_2 T_3;
varexo f E_1 E_2 E_3 U_1 U_2 U_3 V_1 V_2 V_3 W_1 W_2 W_3;
parameters P_1_1 P_1_2 P_1_3 P_2_1 P_2_2 P_2_3 P_3_1 P_3_2 P_3_3;
varexo_det D_1 D_2 D_3;

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
  T_1 = D_1 + W_1;
  T_2 = D_2 + W_2;
  T_3 = D_3 + W_3;
  r = 0.9*r(+1)+f;
end;

shocks;
  var E_1 = 0.1;
  var E_2 = 0.2;
  var E_3 = 0.3;
  var E_1, E_2 = 0.08;

  var U_1 = 0.4;
  var U_2; stderr 0.5;
  var U_3; stderr 0.6;
  corr U_1, U_2 = 0.09;

  var V_1; stderr 1.1;
  var V_2; stderr 1.2;
  var V_3; stderr 1.3;

  var W_1 = 1;
  var W_2 = 5;
  var W_3 = 9;
  var W_1, W_2 = 0.2;
  var W_1, W_3 = 0.3;
  var W_2, W_3 = 0.6;

  var f = 3;
  var f, E_2 = 0.5;
  corr U_1, f = 0.7;
end;

stoch_simul(order=1, nodecomposition, irf=0);

L = load('stochastic_matrix_results.mat');
if max(max(abs(L.oo_.dr.ghu - oo_.dr.ghu))) > 1e-12
  error('Failure in matrix expansion')
end
if max(max(abs(L.oo_.dr.ghx - oo_.dr.ghx))) > 1e-12
  error('Failure in matrix expansion')
end
if max(max(abs(L.oo_.var - oo_.var))) > 1e-12
  error('Failure in matrix expansion')
end
