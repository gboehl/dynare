@#include "fs2000_common.inc"

filter_initial_state;
P(0)=2.24;
k(0)=5.8;
y(0)=0.58;
m(0)=1.01;
end;

varobs gp_obs gy_obs;

estimation(order=1, datafile='../fs2000/fsdat_simul', mode_compute=0,nobs=192, loglinear, smoother, smoothed_state_uncertainty) m P c e W R k d n l gy_obs gp_obs y dA;
estimation(order=1, datafile='../fs2000/fsdat_simul', mode_compute=0,nobs=192, loglinear, smoother,kalman_algo=1, smoothed_state_uncertainty) m P c e W R k d n l gy_obs gp_obs y dA;
estimation(order=1, datafile='../fs2000/fsdat_simul', mode_compute=0,nobs=192, loglinear, smoother,kalman_algo=2, smoothed_state_uncertainty) m P c e W R k d n l gy_obs gp_obs y dA;
estimation(order=1, datafile='../fs2000/fsdat_simul', mode_compute=0,nobs=192, loglinear, smoother,kalman_algo=3, smoothed_state_uncertainty) m P c e W R k d n l gy_obs gp_obs y dA;
estimation(order=1, datafile='../fs2000/fsdat_simul', mode_compute=0,nobs=192, loglinear, smoother,kalman_algo=4, smoothed_state_uncertainty) m P c e W R k d n l gy_obs gp_obs y dA;