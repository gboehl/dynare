// Do not forget to update deterministic_scalar.mod when this file is modified

var(rows=3) X, Y, Z;
var t;
varexo(rows=3) E, U, V;
parameters(rows=3, cols=3) P;

P = [.1, .2, -.3; -.1, .5, .1; .2, -.3, .9];

model;
  X = P*X(-1)+E;
  transpose(Y) = transpose(Y(-1))*transpose(P)+transpose(U);
  Z = V;
  t = Y[1](-1) + Y[2](+1);
end;

shocks;
  var E;
  periods 1 2;
  values [0.5; 0.6; 0.7] [0.2; 0.3; 0.4];

  var U[1];
  periods 3;
  values 0.8;

  var U[2];
  periods 3;
  values 0.7;

  var V[:];
  periods 1;
  values 0.1;
end;

perfect_foresight_setup(periods=50);
perfect_foresight_solver;

write_latex_dynamic_model;
