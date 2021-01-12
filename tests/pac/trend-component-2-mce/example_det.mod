// --+ options: json=compute +--

var x1 x2 x1bar x2bar z y1 y2;

varexo ex1 ex2 ex1bar ex2bar ez ey1 ey2;

parameters a_x1_0 a_x1_1 a_x1_2 a_x1_x2_1 a_x1_x2_2
	   a_x2_0 a_x2_1 a_x2_2 a_x2_x1_1 a_x2_x1_2
	   e_c_m c_z_1 c_z_2 gamma beta ;

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

gamma =  .7;    

trend_component_model(model_name=toto, eqtags=['eq:x1', 'eq:x2', 'eq:x1bar', 'eq:x2bar'], targets=['eq:x1bar', 'eq:x2bar']);

pac_model(discount=beta, model_name=pacman);

model;

[name='eq:exo:1']
diff(y1) = .7*diff(y1(-1)) - .3*diff(y1(-2)) + ey1;

[name='eq:exo:2']
diff(y2) = .5*diff(y2(-1)) - .2*diff(y2(-3)) + ey2;

[name='eq:x1']
diff(x1) = a_x1_0*(x1(-1)-x1bar(-1)) + a_x1_1*diff(x1(-1)) + a_x1_2*diff(x1(-2)) + a_x1_x2_1*diff(x2(-1)) + a_x1_x2_2*diff(x2(-2)) + ex1;

[name='eq:x2']
diff(x2) = a_x2_0*(x2(-1)-x2bar(-1)) + a_x2_1*diff(x1(-1)) + a_x2_2*diff(x1(-2)) + a_x2_x1_1*diff(x2(-1)) + a_x2_x1_2*diff(x2(-2)) + ex2;     

[name='eq:x1bar']
x1bar = x1bar(-1) + ex1bar;

[name='eq:x2bar']
x2bar = x2bar(-1) + ex2bar;

[name='eq:pac']
diff(z) = gamma*(e_c_m*(x1(-1)-z(-1)) + c_z_1*diff(z(-1)) + c_z_2*diff(z(-2)) + pac_expectation(pacman)) + (1-gamma)*(.5*diff(y1)-.7*diff(y2)) + ez;

end;

// Initialize the PAC model.
pac.initialize('pacman');

// Update the PAC/MCE parameters (α in M_.params).
pac.mce.parameters('pacman');

// Setup a scenario for the shocks.
shocks;
   var ex1;
   periods 1;
   values .2;
end;

perfect_foresight_setup(periods=50);
perfect_foresight_solver;