@#include "fs2000.common.inc"

options_.solve_tolf = 1e-12;
estimation(mode_compute=8,silent_optimizer,order=1, datafile='../fs2000/fsdat_simul', nobs=192, mh_replic=0,optim=(
'MaxIter',5000,
'TolFun',1e-4,
'TolX',1e-4,
'MaxFunEvals',5000,
'MaxFunEvalFactor',500,
'InitialSimplexSize',0.05
));
