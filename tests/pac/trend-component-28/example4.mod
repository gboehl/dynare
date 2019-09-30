// --+ options: json=compute, stochastic +--

PLOTS = true;

var x1 x2 x3 x1bar x2bar z y x u v s;

varexo ex1
       ex2
       ex1bar (used='estimationonly')
       ex2bar (used='estimationonly')
       ez
       ey
       ex
       eu
       ev
       es;

parameters
       rho_1 rho_2 rho_3 rho_4
       a_x1_0 a_x1_1 a_x1_2 a_x1_x2_1 a_x1_x2_2
	   a_x2_0 a_x2_1 a_x2_2 a_x2_x1_1 a_x2_x1_2
	   e_c_m c_z_1 c_z_2 c_z_dx2 c_z_u c_z_dv c_z_s cx cy beta
       lambda
       px3;

rho_1 =  .9;
rho_2 = -.2;
rho_3 =  .4;
rho_4 = -.3;


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
c_z_dx2 = .3;
c_z_u = .3;
c_z_dv = .4;
c_z_s  = -.2; 
cx = 1.0;
cy = 1.0;


lambda = 0.5; // Share of optimizing agents.

px3 = -.1;

trend_component_model(model_name=toto, eqtags=['eq:x1', 'eq:x2', 'eq:x1bar', 'eq:x2bar'], targets=['eq:x1bar', 'eq:x2bar']);

pac_model(auxiliary_model_name=toto, discount=beta, model_name=pacman);

model;

[name='eq:u']
s = .3*s(-1) - .1*s(-2) + es;

[name='eq:diff(v)']
diff(v) = .5*diff(v(-1)) + ev;

[name='eq:u']
u = .5*u(-1) - .2*u(-2) + eu;

[name='eq:y']
y = rho_1*y(-1) + rho_2*y(-2) + ey;

[name='eq:x']
x = rho_3*x(-1) + rho_4*x(-2) + ex;

[name='eq:x1']
diff(x1) = a_x1_0*(x1(-1)-x1bar(-1)) + a_x1_1*diff(x1(-1)) + a_x1_2*diff(x1(-2)) + a_x1_x2_1*diff(x2(-1)) + a_x1_x2_2*diff(x2(-2)) + ex1;

[name='eq:x2']
diff(x2) = a_x2_0*(x2(-1)-x2bar(-1)) + a_x2_1*diff(x1(-1)) + a_x2_2*diff(x1(-2)) + a_x2_x1_1*diff(x2(-1)) + a_x2_x1_2*diff(x2(-2)) + ex2;

[name='eq:x3']
diff(log(x3)) = x ;

[name='eq:x1bar']
x1bar = x1bar(-1) + ex1bar;

[name='eq:x2bar']
x2bar = x2bar(-1) + ex2bar;

[name='zpac']
diff(z) = lambda*(e_c_m*(x1(-1)-z(-1)) + c_z_1*diff(z(-1))  + c_z_2*diff(z(-2)) + pac_expectation(pacman) + c_z_s*s + c_z_dv*diff(v) ) + (1-lambda)*( cy*y + cx*x) + c_z_u*u + c_z_dx2*diff(x2) + ez;

end;

shocks;
  var ex1 = 1.0;
  var ex2 = 1.0;
  var ex1bar = 1.0;
  var ex2bar = 1.0;
  var ez = 1.0;
  var ey = 0.1;
  var ex = 0.01;
  var eu = 0.05;
  var ev = 0.05;
  var es = 0.07;
end;

// Initialize the PAC model (build the Companion VAR representation for the auxiliary model).
pac.initialize('pacman');

// Update the parameters of the PAC expectation model (h0 and h1 vectors).
pac.update.expectation('pacman');

// Set initial conditions to zero for non logged variables, and one for logged variables
init = ones(10, M_.endo_nbr+M_.exo_nbr);
initialconditions = dseries(init, 2000Q1, vertcat(M_.endo_names,M_.exo_names));

// Simulate the model for 500 periods
options_.solve_tolf = 1e-9;
TrueData = simul_backward_model(initialconditions, 5000);
TrueData_ = copy(TrueData);

TrueData = equation.evaluate(TrueData, 'eq:x3', 2004Q1);

if PLOTS

    figure(1)
    plot(TrueData_.x.data(12:end), TrueData_.x3.data(12:end)-TrueData.x3.data(12:end), '.k', 'linewidth', 2);

    figure(2)
    plot(TrueData_.x.data(12:end), '-k', 'linewidth', 2);

    figure(3)
    plot([TrueData_.x3.data(12:end)-TrueData.x3.data(12:end)], '-k', 'linewidth', 2);

end

fprintf('Max. abs. error is %s.\n', num2str(max(abs(TrueData.x3.data(12:end)-TrueData_.x3.data(12:end))), 16));

if max(abs(TrueData.x3.data(12:end)-TrueData_.x3.data(12:end)))>1e-5
   error('equation.evaluate() returned wrong values.')
end
