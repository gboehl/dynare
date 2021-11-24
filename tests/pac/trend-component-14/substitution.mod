// --+ options: transform_unary_ops, json=compute, stochastic +--

var x1 x2 x1bar x2bar z y;

varexo ex1 ex2 ex1bar ex2bar ez ey;

parameters
       rho_1 rho_2
       a_x1_0 a_x1_1 a_x1_2 a_x1_x2_1 a_x1_x2_2
	   a_x2_0 a_x2_1 a_x2_2 a_x2_x1_1 a_x2_x1_2
	   e_c_m c_z_1 c_z_2 beta;

rho_1 = .9;
rho_2 = -.2;

a_x1_0 =  -.9;
a_x1_1 =  .4;
a_x1_2 =  .3;
a_x1_x2_1 = .1;
a_x1_x2_2 = .2;

a_x2_0 =  -.9;
a_x2_1 =   .2;
a_x2_2 =  -.1;
a_x2_x1_1 = -.1;
a_x2_x1_2 = .2;

beta  =  .2;
e_c_m =  .5;
c_z_1 =  .2;
c_z_2 = -.1;

@#include "example1/model/pac-expectations/pacman-parameters.inc"

model;

[name='eq:y']
y = rho_1*y(-1) + rho_2*y(-2) + ey;

[name='eq:x1']
diff(x1) = a_x1_0*(x1(-1)-x1bar(-1)) + a_x1_1*diff(x1(-1)) + a_x1_2*diff(x1(-2)) + a_x1_x2_1*diff(x2(-1)) + a_x1_x2_2*diff(x2(-2)) + ex1;

[name='eq:x2']
diff(x2) = a_x2_0*(x2(-1)-x2bar(-1)) + a_x2_1*diff(x1(-1)) + a_x2_2*diff(x1(-2)) + a_x2_x1_1*diff(x2(-1)) + a_x2_x1_2*diff(x2(-2)) + ex2;

[name='eq:x1bar']
x1bar = x1bar(-1) + ex1bar;

[name='eq:x2bar']
x2bar = x2bar(-1) + ex2bar;

[name='zpac']
diff(z) = e_c_m*(x1(-1)-z(-1)) + c_z_1*diff(z(-1))  + c_z_2*diff(z(-2)) +
@#include "example1/model/pac-expectations/pacman-expression.inc"
+ ez;

end;

shocks;
    var ex1 = 1.0;
    var ex2 = 1.0;
    var ex1bar = 1.0;
    var ex2bar = 1.0;
    var ez = 1.0;
    var ey = 0.1;
end;

// Set initial conditions to zero. Please use more sensible values if any...
initialconditions = dseries(zeros(10, M_.endo_nbr+M_.exo_nbr), 2000Q1, vertcat(M_.endo_names,M_.exo_names));

// Simulate the model for 500 periods
set_dynare_seed('default');
TrueData = simul_backward_model(initialconditions, 300);

verbatim;
  set_dynare_seed('default');
  y = zeros(M_.endo_nbr,1);
  y(1:M_.orig_endo_nbr) = rand(M_.orig_endo_nbr, 1);
  x = randn(M_.exo_nbr,1);
  y = substitution.set_auxiliary_variables(y, x, M_.params);
  y = [y(find(M_.lead_lag_incidence(1,:))); y];
  example1 = load('example1.mat');
  [residual, g1] = substitution.dynamic(y, x', M_.params, oo_.steady_state, 1);
end;

if max(abs(example1.TrueData.data(:)-TrueData.data(:)))>1e-9
   error('Simulations do not match.')
end

if ~isequal(length(residual), length(example1.residual)) || max(abs(example1.residual-residual))>1e-8
  warning('Residuals do not match!')
end

if ~isequal(length(g1(:)), length(example1.g1(:))) || max(abs(example1.g1(:)-g1(:)))>1e-8
  warning('Jacobian matrices do not match!')
end

delete('example1.mat');
