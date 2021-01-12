var y x z;

varexo ex ey ez;

parameters a_y_1 a_y_2 b_y_1 b_y_2 b_x_1 b_x_2 c_z_1 c_z_2 c_z_3;

a_y_1 =  .2;
a_y_2 =  .3;
b_y_1 =  .1;
b_y_2 =  .4;
b_x_1 = -.1;
b_x_2 = -.2; 
c_z_1 =  .9;
c_z_2 =  .1;
c_z_3 = -.8;

var_model(model_name=toto, eqtags=['eq:x', 'eq:y', 'eq:z']);

model;

[name='eq:y']
y = a_y_1*y(-1) + a_y_2*diff(x(-1)) + b_y_1*y(-2) + b_y_2*diff(x(-2)) + ey ;

[name='eq:x']
diff(x) = b_x_1*y(-2) + b_x_2*diff(x(-1)) + ex ;

[name='eq:z']
diff(z) = c_z_1*(x(-1)-z(-1)) + c_z_2*diff(z(-1)) + c_z_3*diff(z(-2)) + ez ;

end;

get_ar_matrices('toto');