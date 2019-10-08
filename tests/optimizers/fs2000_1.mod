@#include "fs2000.common.inc"

if (isoctave && user_has_octave_forge_package('optim', '1.6')) || (~isoctave && user_has_matlab_license('optimization_toolbox'))
  estimation(mode_compute=1,order=1, datafile='../fs2000/fsdat_simul', nobs=192, mh_replic=0);
end
