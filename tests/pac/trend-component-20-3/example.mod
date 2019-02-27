// --+ options: json=compute, stochastic +--

// Check that the content of pac.[pacmodel].equations.[eqtag].non_optimizing_behaviour.vars is correct.

var x1 x2 x x1bar x2bar y z ;

varexo ex1 ex2 ex1bar ex2bar ez ey ex;

parameters
       rho_1 rho_2 rho_3 rho_4
       a_x1_0 a_x1_1 a_x1_2 a_x1_x2_1 a_x1_x2_2
	   a_x2_0 a_x2_1 a_x2_2 a_x2_x1_1 a_x2_x1_2
	   e_c_m c_z_1 c_z_2 beta
       lambda;

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

lambda = 0.5; // Share of optimizing agents.

trend_component_model(model_name=toto, eqtags=['eq:x1', 'eq:x2', 'eq:x1bar', 'eq:x2bar'], targets=['eq:x1bar', 'eq:x2bar']);

pac_model(auxiliary_model_name=toto, discount=beta, model_name=pacman);

model;

[name='eq:x1']
diff(x1) = a_x1_0*(x1(-1)-x1bar(-1)) + a_x1_1*diff(x1(-1)) + a_x1_2*diff(x1(-2)) + a_x1_x2_1*diff(x2(-1)) + a_x1_x2_2*diff(x2(-2)) + ex1;     

[name='eq:y']
y = rho_1*y(-1) + rho_2*y(-2) + ey;

[name='eq:x2']
diff(x2) = a_x2_0*(x2(-1)-x2bar(-1)) + a_x2_1*diff(x1(-1)) + a_x2_2*diff(x1(-2)) + a_x2_x1_1*diff(x2(-1)) + a_x2_x1_2*diff(x2(-2)) + ex2;     

[name='eq:x1bar']
x1bar = x1bar(-1) + ex1bar;

[name='eq:x']
x = rho_3*x(-1) + rho_4*x(-2) + ex;

[name='eq:x2bar']
x2bar = x2bar(-1) + ex2bar;

[name='zpac']
diff(z) = lambda*(e_c_m*(x1(-1)-z(-1)) + c_z_1*diff(z(-1))  + c_z_2*diff(z(-2)) + pac_expectation(pacman)) + (1-lambda)*( y + x) + ez;

end;

if ~isequal(M_.endo_names(M_.pac.pacman.equations.eq0.non_optimizing_behaviour.vars), {'y'; 'x'})
   error('PAC non_optimizing_behaviour.vars field is wrong.')
end