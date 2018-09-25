// --+ options: json=compute, stochastic +--

var y x z;

varexo ex ey ez;

parameters a_y_1 a_y_2 b_y_1 b_y_2 b_x_1 b_x_2 ; // VAR parameters

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

var_model(model_name=toto, eqtags=['eq:x', 'eq:y']);

pac_model(auxiliary_model_name=toto, discount=beta, model_name=pacman);

model;

[name='eq:y']
y = a_y_1*y(-1) + a_y_2*diff(x(-1)) + b_y_1*y(-2) + b_y_2*diff(x(-2)) + ey ;

[name='eq:x', data_type='nonstationary']
diff(x) = b_x_1*y(-2) + b_x_2*diff(x(-1)) + ex ;

[name='eq:pac']
diff(z) = e_c_m*(x(-1)-z(-1)) + c_z_1*diff(z(-1))  + c_z_2*diff(z(-2)) + pac_expectation(pacman) + ez;

end;

shocks;
    var ex = 1.0;
    var ey = 1.0;
    var ez = 1.0;
end;

[~, b0, ~] = get_companion_matrix_legacy('toto');

[~, b1, ~] = get_companion_matrix('toto');

if any(abs(b0(:)-b1(:))>1e-9)
   error('get_companion_matrix and get_comapnion_matrix_legacy do not return the same AR matrices.')
end