/*
 * Copyright © 2010-2023 Dynare Team
 *
 * This file is part of Dynare.
 *
 * Dynare is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * Dynare is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with Dynare.  If not, see <https://www.gnu.org/licenses/>.
 */

#ifndef _MEXMOD
#define _MEXMOD

#if defined(MATLAB_MEX_FILE) || defined(OCTAVE_MEX_FILE)

# include <dynmex.h>
# include <dynblas.h>
# include <dynlapack.h>

# define dw_malloc mxMalloc
# define dw_calloc mxCalloc
# define dw_realloc mxRealloc
# define dw_free mxFree
# define dw_exit msExit

/* Handle Ctrl-C in Matlab/Octave */
# ifdef MATLAB_MEX_FILE
extern bool utIsInterruptPending();
# else
#  include <octave/quit.h>
# endif

// NB: C23 has the [[noreturn]] attribute, so this #ifdef can be removed when
// we upgrade
#ifdef __cplusplus
[[noreturn]]
#else
_Noreturn
#endif
void msExit(int status);

extern int constant_seed;

#endif
#endif
