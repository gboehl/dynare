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

options_.simul.maxit = 100;
options_.ep.verbosity = 0;
options_.ep.stochastic.status = 0;
options_.ep.order = 0;
options_.ep.nnodes = 0;
options_.console_mode = 0;

// Extended path simulation
ts = extended_path([], 10, [], options_, M_, oo_);

// Stochastic extended path simulation
options_.ep.stochastic.status = 1;
options_.ep.stochastic.IntegrationAlgorithm='Tensor-Gaussian-Quadrature';
options_.ep.order = 1;
options_.ep.nnodes = 3;
sts = extended_path([], 10, [], options_, M_, oo_);

// The generated paths should be identical (because the model is linear)
if max(max(abs(ts.data-sts.data))) > 1e-12
   error('extended path algorithm fails in ./tests/ep/linearmodel.mod')
end
