//Toy model of Alichi et al. (2015): Avoiding Dark Corners: A Robust Monetary Policy Framework for the United States
// checks whether lmmcp works with ramsey_model

var i y pi;
varexo e_y e_pi;
parameters beta1 beta2 beta3 lambda1 lambda2 pi_bar;

beta1 = 0.6;
beta2 = 0.25;
beta3 = -0.2;
lambda1 = 0.7;
lambda2 = 0.1;
pi_bar = 2.0;

model;
y = beta1*y(-1) + beta2*y(+1) + beta3*(i-pi(+1)) + e_y ⟂ i > 0;
pi = lambda1*pi(+1) + (1-lambda1)*pi(-1) + lambda2*y + e_pi;
end;

planner_objective (pi-pi_bar)^2 + y^2;
ramsey_model(planner_discount=1.0);
histval;
y(0) = -2.0;
pi(0) = 1.0;
end;

steady;

perfect_foresight_setup(periods=50);
perfect_foresight_solver(lmmcp);

rplot i;
