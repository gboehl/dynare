// --+ options: json=compute, stochastic +--

var y x z;

varexo ex ey ez ;

parameters a_y_1 a_y_2 b_y_1 b_y_2 b_x_1 b_x_2 d_y; // VAR parameters

parameters beta e_c_m c_z_1 c_z_2;               // PAC equation parameters

a_y_1 =  .2;
a_y_2 =  .3;
b_y_1 =  .1;
b_y_2 =  .4;
b_x_1 = -.1;
b_x_2 = -.2;
d_y = .5;

beta  =  .9;
e_c_m =  .1;
c_z_1 =  .7;
c_z_2 = -.3;

var_model(model_name=toto, structural, eqtags=['eq:x', 'eq:y']);

pac_model(auxiliary_model_name=toto, discount=beta, model_name=pacman, growth=diff(log(x(-2))));

model;

  [name='eq:y']
  y = a_y_1*y(-1) + a_y_2*diff(log(x(-1))) + b_y_1*y(-2) + b_y_2*diff(log(x(-2))) + ey ;


  [name='eq:x']
  diff(log(x)) = b_x_1*y(-2) + b_x_2*diff(log(x(-1))) + ex ;

  [name='eq:pac']
  diff(log(z)) = e_c_m*(log(x(-1))-log(z(-1))) + c_z_1*diff(log(z(-1)))  + c_z_2*diff(log(z(-2))) + pac_expectation(pacman) + ez;

end;

// Initialize the PAC model (build the Companion VAR representation for the auxiliary model).
pac.initialize('pacman');

// Update the parameters of the PAC expectation model (h0 and h1 vectors).
pac.update.expectation('pacman');

// Print expanded PAC_EXPECTATION term.
pac.print('pacman', 'eq:pac');
