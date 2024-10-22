// --+ options: json=compute, stochastic +--

var y x z v u w z1 z2 z3;

varexo ex ey ez eu ew;

parameters a_y_1 a_y_2 b_y_1 b_y_2 a_u_1 a_u_2 b_u_1 b_u_2 b_x_1 b_x_2 d_y;

parameters beta e_c_m c_z_1 c_z_2;               // PAC equation parameters

a_y_1 =  .2;
a_y_2 =  .3;
b_y_1 =  .1;
b_y_2 =  .4;
a_u_1 =  .2;
a_u_2 =  .3;
b_u_1 =  .1;
b_u_2 =  .4;
b_x_1 = -.1;
b_x_2 = -.2;
d_y = .5;

beta  =  .9;
e_c_m =  .1;
c_z_1 =  .7;
c_z_2 = -.3;

var_model(model_name=toto, structural, eqtags=['eq:x', 'eq:y']);

pac_model(auxiliary_model_name=toto, discount=beta, model_name=pacman);

pac_target_info(pacman);

  target v;
  auxname_target_nonstationary vns;

  component y;
  auxname pv_y_;
  kind ll;

  component x;
  auxname pv_dx_;
  growth diff(x(-1));
  kind dl;

end;

model;

  [name='eq:u']
  u = a_u_1*u(-1) + a_u_2*diff(x(-1)) + b_u_1*y(-2) + b_u_2*diff(x(-2)) + eu ;

  [name='eq:y']
  y = a_y_1*y(-1) + a_y_2*diff(x(-1)) + b_y_1*y(-2) + b_y_2*diff(x(-2)) + ey ;

  [name='eq:x']
  diff(x) = b_x_1*y(-2) + b_x_2*diff(x(-1)) + ex ;

  [name='eq:v']
  v = x + d_y*y ;

  [name='eq:pac']
  diff(z) = e_c_m*(pac_target_nonstationary(pacman)-z(-1)) + c_z_1*diff(z(-1))  + c_z_2*diff(z(-2)) + pac_expectation(pacman) + ez;

  [name='eq:w']
  diff(w) = -.1*w + ew ;

  [name='eq:z1']
  z1 = z+y-x+u;

  [name='eq:z2']
  z2 = z-y+x-u;

  [name='eq:z3']
  z3 = u-diff(v);

end;

shocks;
  var ex = 1.0;
  var ey = 1.0;
  var ez = 1.0;
  var ew = 1.0;
  var eu = 1.0;
end;

// Initialize the PAC model (build the Companion VAR representation for the auxiliary model).
pac.initialize('pacman');

// Update the parameters of the PAC expectation model (h0 and h1 vectors).
pac.update.expectation('pacman');

// Select a subset of the equations and print the equations, the list of parameters, endogenous
// variables and exogenous variables in .inc files under ./simulation-files folder. Note that
// innovations ex1bar and ex2bar will not appear in the equations.
verbatim;
  cherrypick('test2', 'simulation-files-2-a', {'eq:z1', 'eq:z2', 'eq:z3'}, true);
  cherrypick('test2', 'simulation-files-2-b', {'eq:pac', 'eq:x', 'eq:y'}, true);
  aggregate('toto2.mod', {}, '', 'simulation-files-2-a', 'simulation-files-2-b');
end;
