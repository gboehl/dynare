// --+ options: json=compute +--

var y x z;

varexo ex ey ez ;

parameters a_y_1 a_y_2 b_y_1 b_y_2 b_x_1 b_x_2 g; // VAR parameters

parameters beta e_c_m c_z_1 c_z_2;               // PAC equation parameters

a_y_1 =  .2;
a_y_2 =  .3;
b_y_1 =  .1;
b_y_2 =  .4;
b_x_1 = -.1;
b_x_2 = -.2;

beta  =  .9;
e_c_m =  .1;
c_z_1 =  .7;
c_z_2 = -.3;

g = .02;

var_model(model_name=toto, eqtags=['eq:x', 'eq:y']);



model;

[name='eq:y']
y = a_y_1*y(-1) + a_y_2*diff(x(-1)) + b_y_1*y(-2) + b_y_2*diff(x(-2)) + ey ;

[name='eq:x', data_type='nonstationary']
diff(x) = b_x_1*y(-2) + b_x_2*diff(x(-1)) + ex ;

[name='eq:pac']
diff(z) = e_c_m*(x(-1)-z(-1)) + c_z_1*diff(z(-1))  + c_z_2*diff(z(-2)) + pac_expectation(model_name=pacman, var_model_name=toto, discount=beta, growth=g) + ez;

end;

get_companion_matrix('toto');
