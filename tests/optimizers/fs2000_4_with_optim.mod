@#include "fs2000.common.inc"

estimation(mode_compute=4,silent_optimizer,order=1, datafile='../fs2000/fsdat_simul', nobs=192, mh_replic=0,optim=('NumgradEpsilon',1e-6,'NumgradAlgorithm',3,'MaxIter',100,'InitialInverseHessian','eye(9)*.0001'));