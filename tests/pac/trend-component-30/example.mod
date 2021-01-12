// --+ options: json=compute, stochastic +--

// Test with two PAC equations in the same file

var x1 x2 x1bar x2bar z u;

varexo ex1 ex2 ex1bar ex2bar ez eu;

parameters a_x1_0 a_x1_1 a_x1_2 a_x1_x2_1 a_x1_x2_2
	   a_x2_0 a_x2_1 a_x2_2 a_x2_x1_1 a_x2_x1_2
	   e_c_m_z e_c_m_u c_z_1 c_z_2 c_u_1 c_u_2 gamma_z gamma_u beta ;

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
e_c_m_z =  .5;
e_c_m_u =  .6;
c_z_1 =  .2;
c_z_2 = -.1;
c_u_1 =  .2;
c_u_2 = -.1;
gamma_z =  .7;
gamma_u =  .8;

trend_component_model(model_name=toto, eqtags=['eq:x1', 'eq:x2', 'eq:x1bar', 'eq:x2bar'], targets=['eq:x1bar', 'eq:x2bar']);

pac_model(auxiliary_model_name=toto, discount=beta, model_name=pacman);

pac_model(auxiliary_model_name=toto, discount=beta, model_name=nowhereman);

model;

[name='eq:x1']
diff(x1) = a_x1_0*(x1(-1)-x1bar(-1)) + a_x1_1*diff(x1(-1)) + a_x1_2*diff(x1(-2)) + + a_x1_x2_1*diff(x2(-1)) + a_x1_x2_2*diff(x2(-2)) + ex1;     

[name='eq:x2']
diff(x2) = a_x2_0*(x2(-1)-x2bar(-1)) + a_x2_1*diff(x1(-1)) + a_x2_2*diff(x1(-2)) + + a_x2_x1_1*diff(x2(-1)) + a_x2_x1_2*diff(x2(-2)) + ex2;     

[name='eq:x1bar']
x1bar = x1bar(-1) + ex1bar;

[name='eq:x2bar']
x2bar = x2bar(-1) + ex2bar;

[name='eq:pac:z']
diff(z) = gamma_z*(e_c_m_z*(x1(-1)-z(-1)) + c_z_1*diff(z(-1))  + c_z_2*diff(z(-2)) + pac_expectation(pacman)) + (1-gamma_z)*ez;

[name='eq:pac:u']
diff(u) = gamma_u*(e_c_m_u*(x2(-1)-u(-1)) + c_u_1*diff(u(-1))  + c_u_2*diff(u(-2)) + pac_expectation(nowhereman)) + (1-gamma_u)*eu;


end;

shocks;
    var ex1 = 1.0;
    var ex2 = 1.0;
    var ex1bar = 1.0;
    var ex2bar = 1.0;
    var ez = 1.0;
    var eu = 1.0;
end;
