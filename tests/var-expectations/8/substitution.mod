// --+ options: stochastic,transform_unary_ops,json=compute +--

var foo z x y;
varexo e_x e_y e_z;
parameters a b c d e f beta ;

a =  .9;
b = -.2;
c =  .3;
f =  .8;
d =  .5;
e =  .4;

@#include "example/model/var-expectations/varexp-parameters.inc"

beta = 1/(1+.02);

model;
[ name = 'X' ]
diff(log(x)) = a*diff(log(x(-1))) + b*diff(log(x(-2))) + c*diff(z(-2)) + e_x;
[ name = 'Z' ]
diff(z) = f*diff(z(-1)) + e_z;
[ name = 'Y' ]
log(y) = d*log(y(-2)) + e*diff(z(-1)) + e_y;

foo = .5*foo(-1) +
@#include "example/model/var-expectations/varexp-expression.inc"
;
end;

verbatim;
  set_dynare_seed('default');
  y = zeros(M_.endo_nbr,1);
  y(1:M_.orig_endo_nbr) = rand(M_.orig_endo_nbr, 1);
  x = randn(M_.exo_nbr,1);
  y = substitution.set_auxiliary_variables(y, x, M_.params);
  y = [y(find(M_.lead_lag_incidence(1,:))); y];
  [residual, g1] = substitution.dynamic(y, x', M_.params, oo_.steady_state, 1);
  example = load('example.mat');
end;

if max(abs(example.residual-residual))>1e-8
  error('Residuals do not match!')
end

if max(max(abs(example.g1-g1)))>1e-8
  error('Jacobian matrices do not match!')
end

delete('example.mat')