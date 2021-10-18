var k z c U;

varexo e;

parameters beta gamma alpha delta rho s;

beta = 0.987;
gamma = 1;
delta = 0.012;
alpha = 0.4;
rho = 0.95;
s = 0.07; % increased by a factor of 10 from 0.007 to increase risk correction

model;
k=exp(z)*k(-1)^(alpha)-c+(1-delta)*k(-1);
z=rho*z(-1)+s*e;
U=ln(c);
end;

steady_state_model;
z = 0;
U=ln(c);
end;

planner_objective ln(c);
ramsey_model(instruments=(k,c), planner_discount=beta);

initval;
   k = ((1/beta-(1-delta))/alpha)^(1/(alpha-1));
   c = k^alpha-delta*k;
end;

steady;
resid;
