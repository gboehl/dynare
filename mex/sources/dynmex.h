/*
 * Copyright © 2009-2023 Dynare Team
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

#ifndef _DYNMEX_H
#define _DYNMEX_H

#if !defined(MATLAB_MEX_FILE) && !defined(OCTAVE_MEX_FILE)
# error You must define either MATLAB_MEX_FILE or OCTAVE_MEX_FILE
#endif

#include <mex.h>

#if defined(MATLAB_MEX_FILE) && MATLAB_VERSION < 0x0805
# define mxIsScalar(x) (mxGetM(x) == 1 && mxGetN(x) == 1)
#endif

/* The int64_T and uint64_T type are broken under MinGW for MATLAB < R2015b
   (they actually alias long integer types, which are 32-bit) */
#if defined(MATLAB_MEX_FILE) && defined(__MINGW64__) && MATLAB_VERSION < 0x0806
# define int64_T long long
# define uint64_T unsigned long long
#endif
/* NB: the following #ifdef can be removed when we upgrade to C23, since the
   latter has static_assert as a keyword */
#ifdef __cplusplus
static_assert(sizeof(int64_T) == 8, "The int64_T type is buggy");
static_assert(sizeof(uint64_T) == 8, "The uint64_T type is buggy");
#else
_Static_assert(sizeof(int64_T) == 8, "The int64_T type is buggy");
_Static_assert(sizeof(uint64_T) == 8, "The uint64_T type is buggy");
#endif

#endif
