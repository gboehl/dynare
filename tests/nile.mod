var y mu e;
varobs y;
varexo ee emu;
parameters Ve, Vmu, delta, S1, S2;

multinomial( process=A, number_of_regimes=2, probability = [P1]);
multinomial( process=2, number_of_regimes=2, probability = [P1]);
// probability can be a single object representing a vector with a Dirlichet prior (see below) or a vector of calibrated values
//add bernoulli instruction (similar to markov-switching for independent draws) and new options to create bernoulli process
// I don' know if probabilities would be a better name than probability (vector  being implicit)

model;
y = mu + e;
e = ((S1 - 1)*sqrt(delta) + (2 - S1))*sqrt(Ve)*ee;
mu = mu(-1) + (S2-1)*sqrt(Vmu)*emu;
end;

shocks;
var ee; stderr 1;
var emu; stderr 1;
end;

S1.calibration(process=A, regime=1) = 1;
S1.calibration(process=A, regime=2) = 2;
S2.calibration(process=2, regime=1) = 1;
S2.calibration(process=2, regime=2) = 2;

//P2.prior(shape=dirichlet,mean=3,variance=3,params=[16,2]);
P1.prior(shape=dirichlet,mean=3,variance=3,params=[18,2]);

// we need to add Dirichlet shape and params field the number of params must be equal to number of regimes in multinomial declaration
Ve.prior(shape=inv_gamma,mean=6e4,stdev=6,interval=[0,5e4]);
Vmu.prior(shape=inv_gamma,mean=6e4,stdev=6,interval=[0,5e4]);
delta.prior(shape=beta,mean=2,stdev=4,interval=[1,20]);
// I'm not sure what are the parameters used in the paper for the IG and beta distribution

data(file='niledata.csv');
dmm(drop=100,seed=0,thinning_factor=1,replic=500,max_order_of_integration=1,num_nonstationary=1);

// calibrated parameter by regime
//phi1.calibration(regime=1) = 0.2;
//phi1.prior(regime=2,.....);
