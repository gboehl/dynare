/*
 * Copyright © 2007-2024 Dynare Team
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

#ifndef ERROR_HANDLING_HH
#define ERROR_HANDLING_HH

#include <cmath>
#include <sstream>
#include <string>

using namespace std;

struct GeneralException
{
  const string message;
};

struct FloatingPointException : public GeneralException
{
  const string location;
  FloatingPointException(const string& details, string location_arg) :
      GeneralException {"Floating point error: " + details}, location {move(location_arg)}
  {
  }
};

struct UnaryOpException : public FloatingPointException
{
  UnaryOpException(const string& op, double value, string location_arg) :
      FloatingPointException {[=] {
                                // We don’t use std::to_string(), because it uses fixed formatting
                                ostringstream s;
                                s << op << "(X) with X=" << defaultfloat << value;
                                return s.str();
                              }(),
                              move(location_arg)}
  {
  }
};

struct DivideException : public FloatingPointException
{
  DivideException(double v1, double v2, string location_arg) :
      FloatingPointException {[=] {
                                // We don’t use std::to_string(), because it uses fixed formatting
                                ostringstream s;
                                s << "a/X with a=" << defaultfloat << v1 << " and X= " << v2;
                                return s.str();
                              }(),
                              move(location_arg)}
  {
  }
};

struct PowException : public FloatingPointException
{
  PowException(double base, double exponent, string location_arg) :
      FloatingPointException {[=] {
                                // We don’t use std::to_string(), because it uses fixed formatting
                                ostringstream s;
                                s << "X^a with X=" << defaultfloat << base;
                                if (fabs(base) <= 1e-10)
                                  s << " and a=" << exponent;
                                return s.str();
                              }(),
                              move(location_arg)}
  {
  }
};

struct UserException : public GeneralException
{
  UserException() : GeneralException {"User break"}
  {
  }
};

struct FatalException : public GeneralException
{
  FatalException(const string& details) : GeneralException {"Fatal error: " + details}
  {
  }
};

inline void
test_mxMalloc(void* z, int line, const string& file, const string& func, int amount)
{
  if (!z && amount > 0)
    throw FatalException {"mxMalloc: out of memory " + to_string(amount)
                          + " bytes required at line " + to_string(line) + " in function " + func
                          + " (file " + file};
}

#ifdef MATLAB_MEX_FILE
extern "C" bool utIsInterruptPending();
#else
# include <octave/quit.h>
#endif

#endif
