/*
 * Copyright © 2004-2011 Ondra Kamenik
 * Copyright © 2019 Dynare Team
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

#ifndef SYLV_PARAMS_H
#define SYLV_PARAMS_H

#include <ostream>
#include <string>

#if defined(MATLAB_MEX_FILE) || defined(OCTAVE_MEX_FILE)
# include <dynmex.h>
#endif

enum class status
{
  def,
  changed,
  undef
};

template<class _Type>
struct ParamItem
{
protected:
  using _Self = ParamItem<_Type>;
  status s;
  _Type value;

public:
  ParamItem()
  {
    s = status::undef;
  }
  ParamItem(_Type val)
  {
    value = val;
    s = status::def;
  }
  ParamItem(const _Self& item) = default;
  _Self& operator=(const _Self& item) = default;
  _Self&
  operator=(const _Type& val)
  {
    value = val;
    s = status::changed;
    return *this;
  }
  _Type
  operator*() const
  {
    return value;
  }
  status
  getStatus() const
  {
    return s;
  }
  void
  print(std::ostream& out, const std::string& prefix, const std::string& str) const
  {
    if (s == status::undef)
      return;
    out << prefix << str << "= " << value;
    if (s == status::def)
      out << " <default>";
    out << std::endl;
  }
};

class SylvParams
{
public:
  enum class solve_method
  {
    iter,
    recurse
  };

protected:
  class DoubleParamItem : public ParamItem<double>
  {
  public:
    DoubleParamItem() : ParamItem<double>()
    {
    }
    DoubleParamItem(double val) : ParamItem<double>(val)
    {
    }
    DoubleParamItem(const DoubleParamItem& item) = default;
    DoubleParamItem&
    operator=(const double& val)
    {
      ParamItem<double>::operator=(val);
      return *this;
    }
#if defined(MATLAB_MEX_FILE) || defined(OCTAVE_MEX_FILE)
    mxArray* createMatlabArray() const;
#endif
  };

  class IntParamItem : public ParamItem<int>
  {
  public:
    IntParamItem() : ParamItem<int>()
    {
    }
    IntParamItem(int val) : ParamItem<int>(val)
    {
    }
    IntParamItem(const IntParamItem& item) = default;
    IntParamItem&
    operator=(const int& val)
    {
      ParamItem<int>::operator=(val);
      return *this;
    }
#if defined(MATLAB_MEX_FILE) || defined(OCTAVE_MEX_FILE)
    mxArray* createMatlabArray() const;
#endif
  };

  class BoolParamItem : public ParamItem<bool>
  {
  public:
    BoolParamItem() : ParamItem<bool>()
    {
    }
    BoolParamItem(bool val) : ParamItem<bool>(val)
    {
    }
    BoolParamItem(const BoolParamItem& item) = default;
    BoolParamItem&
    operator=(const bool& val)
    {
      ParamItem<bool>::operator=(val);
      return *this;
    }
#if defined(MATLAB_MEX_FILE) || defined(OCTAVE_MEX_FILE)
    mxArray* createMatlabArray() const;
#endif
  };

  class MethodParamItem : public ParamItem<solve_method>
  {
  public:
    MethodParamItem() : ParamItem<solve_method>()
    {
    }
    MethodParamItem(solve_method val) : ParamItem<solve_method>(val)
    {
    }
    MethodParamItem(const MethodParamItem& item) = default;
    MethodParamItem
    operator=(const solve_method& val)
    {
      ParamItem<solve_method>::operator=(val);
      return *this;
    }
#if defined(MATLAB_MEX_FILE) || defined(OCTAVE_MEX_FILE)
    mxArray* createMatlabArray() const;
#endif
  };

public:
  // input parameters
  MethodParamItem method;          // method of solution: iter/recurse
  DoubleParamItem convergence_tol; // norm for what we consider converged
  IntParamItem max_num_iter;       // max number of iterations
  DoubleParamItem bs_norm;         // Bavely Stewart log₁₀ of norm for diagonalization
  BoolParamItem want_check;        // true => allocate extra space for checks
  // output parameters
  BoolParamItem converged;        // true if converged
  DoubleParamItem iter_last_norm; // norm of the last iteration
  IntParamItem num_iter;          // number of iterations
  DoubleParamItem f_err1;         // norm 1 of diagonalization abs. error C−V·F·V⁻¹
  DoubleParamItem f_errI;         // norm ∞ of diagonalization abs. error C−V·F·V⁻¹
  DoubleParamItem viv_err1;       // norm 1 of error I−V·V⁻¹
  DoubleParamItem viv_errI;       // norm ∞ of error I−V·V⁻¹
  DoubleParamItem ivv_err1;       // norm 1 of error I−V⁻¹·V
  DoubleParamItem ivv_errI;       // norm ∞ of error I−V⁻¹·V
  IntParamItem f_blocks;          // number of diagonal blocks of F
  IntParamItem f_largest;         // size of largest diagonal block in F
  IntParamItem f_zeros;           // number of off diagonal zeros in F
  IntParamItem f_offdiag;         // number of all off diagonal elements in F
  DoubleParamItem rcondA1;        // reciprocal cond 1 number of A
  DoubleParamItem rcondAI;        // reciprocal cond ∞ number of A
  DoubleParamItem eig_min;        // minimum eigenvalue of the solved system
  DoubleParamItem mat_err1;       // rel. matrix 1 norm of A·X−B·X·⊗ⁱC−D
  DoubleParamItem mat_errI;       // rel. matrix ∞ norm of A·X−B·X·⊗ⁱC−D
  DoubleParamItem mat_errF;       // rel. matrix Frob. norm of A·X−B·X·⊗ⁱC−D
  DoubleParamItem vec_err1;       // rel. vector 1 norm of A·X−B·X·⊗ⁱC−D
  DoubleParamItem vec_errI;       // rel. vector ∞ norm of A·X−B·X·⊗ⁱC−D
  DoubleParamItem cpu_time;       // time of the job in CPU seconds

  SylvParams(bool wc = false) :
      method(solve_method::recurse), convergence_tol(1.e-30), max_num_iter(15), bs_norm(1.3),
      want_check(wc)
  {
  }
  SylvParams(const SylvParams& p) = default;
  SylvParams& operator=(const SylvParams& p) = default;
  ~SylvParams() = default;
  void print(const std::string& prefix) const;
  void print(std::ostream& fdesc, const std::string& prefix) const;
  void setArrayNames(int& num, const char** names) const;
#if defined(MATLAB_MEX_FILE) || defined(OCTAVE_MEX_FILE)
  mxArray* createStructArray() const;
#endif
private:
  void copy(const SylvParams& p);
};

inline std::ostream&
operator<<(std::ostream& out, SylvParams::solve_method m)
{
  switch (m)
    {
    case SylvParams::solve_method::iter:
      out << "iterative";
      break;
    case SylvParams::solve_method::recurse:
      out << "recurse (a.k.a. triangular)";
      break;
    }
  return out;
}

#endif /* SYLV_PARAMS_H */
