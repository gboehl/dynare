/*
 * Copyright (C) 2014 Dynare Team
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
 * along with Dynare.  If not, see <http://www.gnu.org/licenses/>.
 */

#include <dynmex.h>
#include <algorithm>
#include <string.h>
using namespace std;

#if defined(__CYGWIN32__) || defined(_WIN32)
extern "C" void dmmmain(char *nmlfile, size_t len);
#else
extern "C" void dmmmain_(char *nmlfile, size_t len);
#endif

void
mexFunction(int nlhs, mxArray *plhs[],
            int nrhs, const mxArray *prhs[])
{
  char nmlfilec[200];
  char nmlfilef[200];

  /*
   * Check args
   */
  if (nrhs != 1 || !mxIsChar(prhs[0]) || nlhs != 0)
    DYN_MEX_FUNC_ERR_MSG_TXT("Error in DMM MEX file: this function takes 1 string input argument and returns no output arguments.");

  /*
   * Get nml file name
   */
  if (mxGetString(prhs[0], nmlfilec, mxGetN(prhs[0])+1))
    DYN_MEX_FUNC_ERR_MSG_TXT("Error in DMM MEX file: error using mxGetString.\n");

  /*
   * Convert to Fortran
   */
  copy(nmlfilec, nmlfilec+strlen(nmlfilec), nmlfilef);
  fill(nmlfilef + strlen(nmlfilec), nmlfilef + 200, ' ');

  /*
   * Call top_level function (formerly main)
   */
  try
    {
#if defined(__CYGWIN32__) || defined(_WIN32)
      dmmmain(nmlfilef, 200);
#else
      dmmmain_(nmlfilef, 200);
#endif
    }
  catch (const char *str)
    {
      DYN_MEX_FUNC_ERR_MSG_TXT(str);
    }
}
