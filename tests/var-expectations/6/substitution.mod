// --+ options: stochastic,transform_unary_ops,json=compute +--

var foo x1 x2 x1bar x2bar;

varexo ex1 ex2 ex1bar ex2bar;

parameters a_x1_0 a_x1_0_ a_x1_1 a_x1_2 a_x1_x2_1 a_x1_x2_2
	   a_x2_0 a_x2_1 a_x2_2 a_x2_x1_1 a_x2_x1_2
       beta ;

a_x1_0 =  -.9;
a_x1_0_ =  -.8;
a_x1_1 =  .4;
a_x1_2 =  .3;
a_x1_x2_1 = .1;
a_x1_x2_2 = .2;


a_x2_0 =  -.9;
a_x2_1 =   .2;
a_x2_2 =  -.1;
a_x2_x1_1 = -.1;
a_x2_x1_2 = .2;

@#include "example/model/var-expectations/varexp-parameters.inc"

beta = 1/(1+.02);

model;

[name='eq:x1', data_type='nonstationary']
diff(x1) = a_x1_0*(x1(-1)-x1bar(-1))+a_x1_0_*(x2(-1)-x2bar(-1)) + a_x1_1*diff(x1(-1)) + a_x1_2*diff(x1(-2)) + + a_x1_x2_1*diff(x2(-1)) + a_x1_x2_2*diff(x2(-2)) + ex1;

[name='eq:x2', data_type='nonstationary']
diff(x2) = a_x2_0*(x2(-1)-x2bar(-1)) + a_x2_1*diff(x1(-1)) + a_x2_2*diff(x1(-2)) + + a_x2_x1_1*diff(x2(-1)) + a_x2_x1_2*diff(x2(-2)) + ex2;

[name='eq:x1bar', data_type='nonstationary']
x1bar = x1bar(-1) + ex1bar;

[name='eq:x2bar', data_type='nonstationary']
x2bar = x2bar(-1) + ex2bar;

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
  example = load('example.mat');
  [residual, g1] = substitution.dynamic(y, x', M_.params, oo_.steady_state, 1);
end;

if max(abs(example.residual-residual))>1e-8
  error('Residuals do not match!')
end

if max(max(abs(example.g1-g1)))>1e-8
  error('Jacobian matrices do not match!')
end

delete('example.mat');