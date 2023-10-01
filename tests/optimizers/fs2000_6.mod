@#include "fs2000.common.inc"

estimation(mode_compute=6,silent_optimizer,order=1, datafile='../fs2000/fsdat_simul', nobs=192, mh_replic=0, optim=('nclimb-mh', 10, 'ncov-mh', 1000, 'nscale-mh', 5000));

// test the mode file generated with mode_compute=6
estimation(order=1,silent_optimizer,datafile='../fs2000/fsdat_simul',nobs=192,loglinear,mode_compute=0,mode_file='fs2000_6/Output/fs2000_6_mode',mh_replic=10,
        posterior_sampler_options=('scale_file','fs2000_6/Output/fs2000_6_optimal_mh_scale_parameter'));
