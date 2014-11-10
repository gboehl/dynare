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

extern "C" void dmmmain_(char *nmlfile, size_t len);

void
mexFunction(int nlhs, mxArray *plhs[],
            int nrhs, const mxArray *prhs[])
{
  char nmlfile[200];

  /*
   * Check args
   */
  if (nrhs != 1 || !mxIsChar(prhs[0]) || nlhs != 0)
    DYN_MEX_FUNC_ERR_MSG_TXT("Error in DMM MEX file: this function takes 1 string input argument and returns no output arguments.");

  /*
   * Get nml file name
   */
  if (mxGetString(prhs[0], nmlfile, mxGetN(prhs[0])+1))
    DYN_MEX_FUNC_ERR_MSG_TXT("Error in DMM MEX file: error using mxGetString.\n");

  /*
   * Call top_level function (formerly main)
   */
  try
    {
      dmmmain_(nmlfile, 200);
    }
  catch (const char *str)
    {
      DYN_MEX_FUNC_ERR_MSG_TXT(str);
    }
}
