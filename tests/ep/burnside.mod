addpath ..
var y x;
varexo e;

parameters beta theta rho xbar;
xbar = 0.0179;
rho =  -0.139;
theta = -1.5;
beta = 0.95;

model(use_dll);
1 = beta*exp(theta*x(+1))*(1+y(+1))/y;
x = (1-rho)*xbar + rho*x(-1)+e;
end;

shocks;
var e; stderr 0.0348;
end;

initval;
x = xbar;
y = beta*exp(theta*xbar)/(1-beta*exp(theta*xbar));
end;

resid;

steady;

if beta*exp(theta*xbar+.5*theta^2*M_.Sigma_e/(1-rho)^2)>1-eps
   disp('The model doesn''t have a solution!')
   return
end

set_dynare_seed('default');
stoch_simul(order=1,irf=0,periods=201);
y_perturbation_1 = oo_.endo_simul(1,:)';

set_dynare_seed('default');
stoch_simul(order=2,irf=0,periods=201);
y_perturbation_2 = oo_.endo_simul(1,:)';

set_dynare_seed('default');
stoch_simul(order=2,pruning,irf=0,periods=201);
y_perturbation_2_pruning = oo_.endo_simul(1,:)';

options_.simul.maxit = 100;

set_dynare_seed('default');
ytrue=exact_solution(M_,oo_, 800);


set_dynare_seed('default');
options_.ep.stochastic.order = 0;
ts = extended_path([], 200, [], options_, M_, oo_);

set_dynare_seed('default');
options_.ep.stochastic.order = 2;
options_.ep.stochastic.IntegrationAlgorithm='Stroud-Cubature-3';//'Unscented'; //'Tensor-Gaussian-Quadrature';
options_.ep.stochastic.quadrature.nodes = 3;
ts1_4 = extended_path([], 200, [], options_, M_, oo_);

set_dynare_seed('default');
options_.ep.stochastic.order = 4;
options_.ep.stochastic.IntegrationAlgorithm='Stroud-Cubature-3';//'Unscented'; //'Tensor-Gaussian-Quadrature';
options_.ep.stochastic.quadrature.nodes = 3;
ts1_8 = extended_path([], 200, [], options_, M_, oo_);

set_dynare_seed('default');
options_.ep.stochastic.order = 6;
options_.ep.stochastic.IntegrationAlgorithm='Stroud-Cubature-3';//'Unscented'; //'Tensor-Gaussian-Quadrature';
options_.ep.stochastic.quadrature.nodes = 3;
ts1_12 = extended_path([], 200, [], options_, M_, oo_);

set_dynare_seed('default');
options_.ep.stochastic.order = 8;
options_.ep.stochastic.IntegrationAlgorithm='Stroud-Cubature-3';//'Unscented'; //'Tensor-Gaussian-Quadrature';
options_.ep.stochastic.quadrature.nodes = 3;
ts1_16 = extended_path([], 200, [], options_, M_, oo_);

disp('True mean and standard deviation')
disp(mean(ytrue(101:end)))
disp(sqrt(var(ytrue(101:end))))

disp('Order 1 perturbation mean and standard deviation')
disp(mean(y_perturbation_1(101:end)))
disp(sqrt(var(y_perturbation_1(101:end))))

disp('Order 2 perturbation mean and standard deviation')
disp(mean(y_perturbation_2(101:end)))
disp(sqrt(var(y_perturbation_2(101:end))))

disp('Order 2 perturbation (with pruning) mean and standard deviation')
disp(mean(y_perturbation_2_pruning(101:end)))
disp(sqrt(var(y_perturbation_2_pruning(101:end))))

disp('Extended path mean and standard deviation')
disp(mean(ts.data(101:end,1)))
disp(sqrt(var(ts.data(101:end,1))))

disp('Stochastic extended path mean and standard deviation')
disp(mean(ts1_4.data(101:end,1)))
disp(sqrt(var(ts1_4.data(101:end,1))))

disp('Accuracy error (order 1 perturbation)')
disp(mean(100*abs(y_perturbation_1-ytrue')./ytrue'));
disp(max(100*abs(y_perturbation_1-ytrue')./ytrue'));

disp('Accuracy error (order 2 perturbation)')
disp(mean(100*abs(y_perturbation_2-ytrue')./ytrue'));
disp(max(100*abs(y_perturbation_2-ytrue')./ytrue'));

disp('Accuracy error (order 2 perturbation with pruning)')
disp(100*mean(abs(y_perturbation_2_pruning-ytrue')./ytrue'));
disp(100*max(abs(y_perturbation_2_pruning-ytrue')./ytrue'));

disp('Accuracy error (extended path)')
disp(mean(100*abs(ts.data(:,1)-ytrue')./ytrue'));
disp(max(100*abs(ts.data(:,1)-ytrue')./ytrue'));

disp('Accuracy error (stochastic extended path (order=2))')
disp(mean(100*abs(ts1_4.data(:,1)-ytrue')./ytrue'));
disp(max(100*abs(ts1_4.data(:,1)-ytrue')./ytrue'));

disp('Accuracy error (stochastic extended path (order=4))')
disp(mean(100*abs(ts1_8.data(:,1)-ytrue')./ytrue'));
disp(max(100*abs(ts1_8.data(:,1)-ytrue')./ytrue'));

disp('Accuracy error (stochastic extended path (order=6))')
disp(mean(100*abs(ts1_12.data(:,1)-ytrue')./ytrue'));
disp(max(100*abs(ts1_12.data(:,1)-ytrue')./ytrue'));

disp('Accuracy error (stochastic extended path (order=8))')
disp(mean(100*abs(ts1_16.data(:,1)-ytrue')./ytrue'));
disp(max(100*abs(ts1_16.data(:,1)-ytrue')./ytrue'));
