%mod-file triggering the sim1_linear.m solver;
%The exogenous arma processes test whether the Jacobian at the
%deterministic steady state is correctly computed
var x
    y
    z;

varexo u
       v;

parameters a1 a2 a3 a4
	   b1 b2 b3
	   c1;

a1 =  .50;
a2 =  .00;
a3 =  .70;
a4 =  .40;
b1 =  .90;
b2 =  .00;
b3 =  .80;
c1 =  .95;


model(linear);
   y = a1*x(-1) + a2*x(+1) + a3*z + a4*y(-1);
   z = b1*z(-1) + b2*z(+1) + b3*x + u;
   x = c1*x(-1) + v +v(-1)+v(+1);
end;

initval;
y=-1;
x=-1;
z=-1;
end;

endval;
y=0;
x=0;
z=0;
end;
steady;
simul(periods=1000,stack_solve_algo=0);
