var x y z;

varexo ex ey ez;

parameters a_y_1 a_y_2 b_y_1 b_y_2 b_x_1 b_x_2 c_z_1 d;

a_y_1 =  .2;
a_y_2 =  .3;
b_y_1 =  .1;
b_y_2 =  .4;
b_x_1 = -.1;
b_x_2 = -.2; 
c_z_1 =  .9;
d     =  .1;

var_model(model_name=toto, eqtags=['eq:x', 'eq:y']);

model;

[name='eq:y']
y = a_y_1*y(-1) + a_y_2*x(-1) + b_y_1*y(-2) + b_y_2*x(-2) + ey ;

z = c_z_1*z(-3) + d*y + ez;

[name='eq:x']
x = b_x_1*y(-2) + b_x_2*x(-2) + ez ;

end;

get_ar_matrices('toto');
