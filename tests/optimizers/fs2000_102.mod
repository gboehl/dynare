@#include "fs2000.common.inc"

if ~isoctave() && exist('simulannealbnd','file')
   estimation(mode_compute=102,silent_optimizer,mode_file='../estimation/fs2000/Output/fs2000_mode',order=1, datafile='../fs2000/fsdat_simul', nobs=192, mh_replic=0, mh_nblocks=2, mh_jscale=0.8);
end
