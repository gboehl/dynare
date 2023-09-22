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

#if defined(MATLAB_MEX_FILE) || defined(OCTAVE_MEX_FILE)

# include <dynmex.h>

# ifdef __cplusplus
extern "C"
{
# endif

  int constant_seed;

  [[noreturn]]
  void
  msExit([[maybe_unused]] int status)
  {
    throw "Error in MS-SBVAR MEX file.\n";
  }

# ifdef __cplusplus
}
# endif
#endif
