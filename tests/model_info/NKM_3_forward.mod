parameters betta sigma phipi phiy omega kappa r_ss;

var c r ui;

varexo pai;

betta = 0.99;         % discount rate patient
sigma = 1.5;          % risk aversion param
omega = 0.75;         % calvo param
phipi = 1.5;
phiy = 0.5;
r_ss = 0;
kappa = (1-omega)*(1-omega*betta)/omega;  

model(linear);
[name='Euler']
c = c(+1) - (1/sigma)*(r - r_ss - pai(+1));
[name='NKPC']
pai = betta*pai(+1) + kappa*c;
[name='Taylor rule']
r = phipi*pai + phiy*c + ui;

end;

initval;
ui = 0;
c = 0;
r = 0;
end;

steady;

shocks;
var pai;
periods 1,2,3,4,5,6,7,8,9,10;
values 0.006,	0.005854795, 0.00550196, 0.0052, 0.0048, 0.0043, 0.0037, 0.0028, 0.0022, 0.0015;
end;
model_info;
model_info(block_static,incidence);

model_info(block_dynamic,incidence);
