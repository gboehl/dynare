// --+ options: json=compute, stochastic +--

var x1 x2 x1bar x2bar z y x u v s;

varexo ex1 ex2 ex1bar ex2bar ez ey ex eu ev es;

parameters
       rho_1 rho_2 rho_3 rho_4
       a_x1_0 a_x1_1 a_x1_2 a_x1_x2_1 a_x1_x2_2
	   a_x2_0 a_x2_1 a_x2_2 a_x2_x1_1 a_x2_x1_2
	   e_c_m c_z_1 c_z_2 c_z_dx2 c_z_u c_z_dv c_z_s cx cy beta
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
c_z_dx2 = .3;
c_z_u = .3;
c_z_dv = .4;
c_z_s  = -.2; 
cx = 1.0;
cy = 1.0;


lambda = 0.5; // Share of optimizing agents.

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
  var ex = 0.1;
  var eu = 0.05;
  var ev = 0.05;
  var es = 0.07;
end;

// Initialize the PAC model (build the Companion VAR representation for the auxiliary model).
pac.initialize('pacman');

// Update the parameters of the PAC expectation model (h0 and h1 vectors).
pac.update.expectation('pacman');

// Exogenous variables with non zero mean:
pac.bgp.set('pacman', 'zpac', 'c_z_s', true);
pac.bgp.set('pacman', 'zpac', 'cx', true);
pac.bgp.set('pacman', 'zpac', 'c_z_dx2', true);

id = find(strcmp('c_z_s', M_.param_names));
id = find(id==M_.pac.pacman.equations.eq0.additive.params);
if ~pac.bgp.get('pacman', 'zpac', 'additive', id)
   error('bgp field in additive is wrong.')
end

id = find(strcmp('c_z_dv', M_.param_names));
id = find(id==M_.pac.pacman.equations.eq0.additive.params);
if pac.bgp.get('pacman', 'zpac', 'additive', id)
   error('bgp field in additive is wrong.')
end

id = find(strcmp('cx', M_.param_names));
id = find(id==M_.pac.pacman.equations.eq0.non_optimizing_behaviour.params);
if ~pac.bgp.get('pacman', 'zpac', 'non_optimizing_behaviour', id)
   error('bgp field in non_optimizing_behaviour is wrong.')
end

id = find(strcmp('cy', M_.param_names));
id = find(id==M_.pac.pacman.equations.eq0.non_optimizing_behaviour.params);
if pac.bgp.get('pacman', 'zpac', 'non_optimizing_behaviour', id)
   error('bgp field in non_optimizing_behaviour is wrong.')
end

id = find(strcmp('c_z_dx2', M_.param_names));
id = find(id==M_.pac.pacman.equations.eq0.optim_additive.params);
if ~pac.bgp.get('pacman', 'zpac', 'optim_additive', id)
   error('bgp field in optim_additive is wrong.')
end

id = find(strcmp('c_z_u', M_.param_names));
id = find(id==M_.pac.pacman.equations.eq0.additive.params);
if pac.bgp.get('pacman', 'zpac', 'optim_additive', id)
   error('bgp field in optim_additive is wrong.')
end