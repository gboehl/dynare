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

#include <cstdlib>
#include <cstring>
#include <cctype>
#include <dynmex.h>

#include "modify_for_mex.h"

int main(int nargs, char **args);

void
mexFunction(int nlhs, [[maybe_unused]] mxArray *plhs[],
            int nrhs, const mxArray *prhs[])
{
  int nargs = 0;
  const char *mainarg = "./a.out";

  /*
   * Check args
   */
  if (nrhs != 1 || !mxIsChar(prhs[0]) || nlhs != 0)
    mexErrMsgTxt("Error in MS-SBVAR MEX file: this function takes 1 string input argument and returns no output argument.");

  /*
   * Allocate memory
   */
  int maxnargs = static_cast<int>(mxGetN(prhs[0])/2+1);
  char *argument = static_cast<char *>(mxCalloc(mxGetN(prhs[0])+1, sizeof(char)));
  char **args = static_cast<char **>(mxCalloc(maxnargs, sizeof(char *)));
  if (!argument || !args)
    mexErrMsgTxt("Error in MS-SBVAR MEX file: could not allocate memory. (1)");

  /*
   * Create argument string from prhs and parse to create args / nargs
   */
  if (!(args[nargs] = static_cast<char *>(mxCalloc(strlen(mainarg)+1, sizeof(char)))))
    mexErrMsgTxt("Error in MS-SBVAR MEX file: could not allocate memory. (2)");

  strcpy(args[nargs++], mainarg);

  if (mxGetString(prhs[0], argument, mxGetN(prhs[0])+1))
    mexErrMsgTxt("Error in MS-SBVAR MEX file: error using mxGetString.\n");

  char *beginarg = argument;
  while (int n = strcspn(beginarg, " "))
    {
      if (!(args[nargs] = static_cast<char *>(mxCalloc(n+1, sizeof(char)))))
        mexErrMsgTxt("Error in MS-SBVAR MEX file: could not allocate memory. (3)");

      strncpy(args[nargs++], beginarg, n);
      beginarg += (isspace(beginarg[n]) || isblank(beginarg[n]) ? ++n : n);
    }
  mxFree(argument);

  /*
   * Call top_level function (formerly main)
   */
  try
    {
      main(nargs, args);
    }
  catch (const char *str)
    {
      mexErrMsgTxt(str);
    }

  /*
   * free memory
   */
  for (int n = 0; n < nargs; n++)
    mxFree(args[n]);
  mxFree(args);
}
