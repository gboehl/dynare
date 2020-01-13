@#include "fs2000.common.inc"

% Skip the test under R2009b, because fminunc fails due to Inf value.
% It remains to be determined in which version it started to work.
if exist('fminunc','file') && (isoctave || ~matlab_ver_less_than('7.10'))
  estimation(mode_compute=3,order=1, datafile='../fs2000/fsdat_simul', nobs=192, mh_replic=0);
end
