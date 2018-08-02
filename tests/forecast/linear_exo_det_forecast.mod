// Regression test for issue #1608 (combining forecast, varexo_det and linear)

var y, pi, i, g, u, k;
varexo e_g e_u e_k;
varexo_det gov;

parameters lambda, pi_target, y_target, phi_pi, phi, rho, rhoout, rhopi, rhoint, sigma1, sigma2, sigma3;

lambda     = 0.3;
pi_target  = 0;
y_target   = 0;
phi_pi      = 1.5;
phi          = 1;
rho          = 0.99;
T            = 50;
rhoout     = 0.8;
rhopi       = 0.5;
rhoint      = 0;
sigma1    = 1;
sigma2    = 1;
sigma3    = 1;

model(linear);
y=y(+1)-phi*(i-pi(+1))+gov+g;
pi=lambda*y+rho*pi(+1)+u;
i=phi_pi*(pi-pi_target)+k;
g=rhoout*g(-1)+e_g;
u=rhopi*u(-1)+e_u;
k=rhoint*k(-1)+e_k;
end;

steady;
check;

shocks;
var e_g;
stderr sigma1;
var e_u;
stderr sigma2;
var e_k;
stderr sigma3;
var gov;
periods 1:9;
values 0.2;
end;

stoch_simul(irf=0);
forecast;
