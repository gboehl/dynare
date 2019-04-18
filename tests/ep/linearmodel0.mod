var y pie r;
varexo e_y e_pie;

parameters delta sigma alpha kappa gamma1 gamma2;

delta =  0.44;
kappa =  0.18;
alpha =  0.48;
sigma = -0.06;

gamma1 = 1.5;
gamma2 = 0.5;

model(use_dll);
y  = delta * y(-1)  + (1-delta)*y(+1)+sigma *(r - pie(+1)) + e_y;
pie  =   alpha * pie(-1) + (1-alpha) * pie(+1) + kappa*y + e_pie;
r = gamma1*pie+gamma2*y;
end;

shocks;
var e_y;
stderr 0.63;
var e_pie;
stderr 0.4;
end;

steady;

stoch_simul(order=1,periods=200, irf=0, nomoments);

oo = oo_;

options_.simul.maxit = 100;
options_.ep.order = 0;

ts = extended_path([], 200, oo_.exo_simul, options_, M_, oo_);

if max(max(abs(transpose(ts.data(2:end,:))-oo.endo_simul))) > 1e-6
   error('extended path algorithm fails in ./tests/ep/linearmodel.mod')
end
