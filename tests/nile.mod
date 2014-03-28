var y e mu;
varexo ee emu;


parameters delta S1 S2 Ve Vmu;

bernoulli( process=1, number_of_regimes=2, parameters = [S1], probability = P1, values{optional}=[1,3]); // first regime = 1, .....
bernoulli( process=2, number_of_regimes=2, parameters = [S2], probability = P2);
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

estimated_params;
probability P1.prior(shape=dirichlet,params=[16,2]);
probability P2.prior(shape=dirichlet,params=[16,2]);
// we need to add Dirichlet shape and params field the number of params must be equal to number of regimes in bernoulli declaration
Ve.prior(shape=inverted_gamma,mean=??,variance=??,bounds=[0,5e4]);
Vmu.prior(shape=inverted_gamma,mean=??,variance=??,bounds=[0,5e4]);
delta.prior(shape=beta,mean=??,variance=??,bounds=[1,20]);
// I'm not sure what are the parameters used in the paper for the IG and beta distribution
end;

dmm(......);
// need to add all DMM options

// calibrated parameter by regime
phi1.calibration(regime=1) = 0.2;
phi1.prior(regime=2,.....);