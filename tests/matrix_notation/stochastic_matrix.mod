// Do not forget to update stochastic_scalar.mod when this file is modified

var(rows=3) X, Y, Z, T;
var r;
varexo(rows=3) E, U, V, W;
varexo f;
parameters(rows=3, cols=3) P;
varexo_det(rows=3) D;

P = [.1, .2, -.3; -.1, .5, .1; .2, -.3, .9];

model;
  X = P*X(-1)+E;
  transpose(Y) = transpose(Y(-1))*transpose(P)+transpose(U);
  Z[1] = V[1];
  Z[2] = V[2];
  Z[3] = V[3];
  T = D + W;
  r = 0.9*r(+1)+f;
end;

shocks;
  var E = [ 0.1, 0.2, 0.3 ];
  var E[1], E[2] = 0.08;

  var U[1] = 0.4;
  var U[2]; stderr 0.5;
  var U[3]; stderr 0.6;
  corr U[1], U[2] = 0.09;

  var V; stderr [ 1.1, 1.2, 1.3 ];

  vcov W = [
      1, 0.2, 0.3;
    0.2, 5,   0.6;
    0.3, 0.6,   9
  ];

  var f = 3;
  var f, E[2] = 0.5;
  corr U[1], f = 0.7;
end;

stoch_simul(order=1, nodecomposition, irf=0);

write_latex_dynamic_model;
