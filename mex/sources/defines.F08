! Under Octave, MX_HAS_INTERLEAVED_COMPLEX is defined through the command-line
! defines (as is done by mkoctfile), while under MATLAB it is defined through
! the present file (because it is defined by mex.h there).

#ifdef MATLAB_MEX_FILE
# define MX_HAS_INTERLEAVED_COMPLEX 1
#endif
