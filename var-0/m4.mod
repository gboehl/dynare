var y x z;

varexo ex ey ez;

parameters a_x_0 a_x_1 a_x_2 b_z_0 b_z_1 b_z_2;

a_x_0 =  .2;
a_x_1 =  .9;
a_x_2 = -.2;
b_z_0 =  .3;
b_z_1 =  .7;
b_z_2 = -.4;

var_model(model_name=toto, eqtags=['eq:x', 'eq:z', 'eq:y']);

model;

[name='eq:y']
y = y(-1) + ey ;

[name='eq:x']
diff(x) = -a_x_0*(x(-1)-y(-1)) + a_x_1*diff(x(-1)) + a_x_2*diff(x(-2)) + ex ;

[name='eq:z']
diff(z) = -b_z_0*(z(-1)-x(-1)) + b_z_1*diff(z(-1)) + b_z_2*diff(z(-2)) + ez ;

end;

get_ar_matrices('toto');