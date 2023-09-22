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

#ifndef _MODIFY_FOR_MEX_H
#define _MODIFY_FOR_MEX_H

#include <dynmex.h>
#include <dynblas.h>
#include <dynlapack.h>

#define dw_malloc mxMalloc
#define dw_calloc mxCalloc
#define dw_realloc mxRealloc
#define dw_free mxFree
#define dw_exit msExit

/* Handle Ctrl-C in Matlab/Octave */
#ifdef MATLAB_MEX_FILE
# ifdef __cplusplus
extern "C"
# else
extern
# endif
bool utIsInterruptPending();
#else
# include <octave/quit.h>
#endif

#ifdef __cplusplus
extern "C"
#endif
// NB: C23 has the [[noreturn]] attribute, so this #ifdef can be removed when
// we upgrade
#ifdef __cplusplus
[[noreturn]]
#else
_Noreturn
#endif
void msExit(int status);

#endif // _MODIFY_FOR_MEX_H
