// --+ options: transform_unary_ops, json=compute, stochastic +--

var y x z;

varexo ex ey ez;

parameters a_y_1 a_y_2 b_y_1 b_y_2 b_x_1 b_x_2 g; // VAR parameters

parameters e_c_m c_z_1 c_z_2;               // PAC equation parameters

a_y_1 =  .2;
a_y_2 =  .3;
b_y_1 =  .1;
b_y_2 =  .4;
b_x_1 = -.1;
b_x_2 = -.2;

g = .1;

e_c_m =  .1;
c_z_1 =  .7;
c_z_2 = -.3;

@#include "example/model/pac-expectations/eq0-pacman-parameters.inc"

model;

[name='eq:y']
y = a_y_1*y(-1) + a_y_2*diff(x(-1)) + b_y_1*y(-2) + b_y_2*diff(x(-2)) + ey ;

[name='eq:x']
diff(x) = b_x_1*y(-2) + b_x_2*diff(x(-1)) + g*(1-b_x_2)  + ex ;

[name='eq:pac']
diff(z) = e_c_m*(x(-1)-z(-1)) + c_z_1*diff(z(-1))  + c_z_2*diff(z(-2)) +
@#include "example/model/pac-expectations/eq0-pacman-growth-neutrality-correction.inc"
+
@#include "example/model/pac-expectations/eq0-pacman-expression.inc"
+ ez;

end;

shocks;
    var ex = 1.0;
    var ey = 1.0;
    var ez = 1.0;
end;

// Set initial conditions to zero. Please use more sensible values if any...
initialconditions = dseries(zeros(10, M_.endo_nbr+M_.exo_nbr), 2000Q1, vertcat(M_.endo_names,M_.exo_names));

// Simulate the model for 20 periods
set_dynare_seed('default');
TrueData = simul_backward_model(initialconditions, 20);

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

if max(abs(example.TrueData.data(:)-TrueData.data(:)))>1e-9
   error('Simulations do not match.')
end

if ~isequal(length(residual), length(example.residual)) || max(abs(example.residual-residual))>1e-8
  warning('Residuals do not match!')
end

if ~isequal(length(g1(:)), length(example.g1(:))) || max(abs(example.g1(:)-g1(:)))>1e-8
  warning('Jacobian matrices do not match!')
end

delete('example.mat');