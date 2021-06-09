/*
 * Copyright © 2005 Ondra Kamenik
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

// Vector function.

/* This file defines interface for functions taking a vector as an input and
   returning a vector (with a different size) as an output. We are also
   introducing a parameter signalling; it is a boolean vector which tracks
   parameters which were changed from the previous call. The VectorFunction
   implementation can exploit this information and evaluate the function more
   efficiently. The information can be completely ignored.

   From the signalling reason, and from other reasons, the function evaluation
   is not const. */

#ifndef VECTOR_FUNCTION_H
#define VECTOR_FUNCTION_H

#include "Vector.hh"
#include "GeneralMatrix.hh"

#include <vector>
#include <memory>

/* This is a simple class representing a vector of booleans. The items night be
   retrieved or changed, or can be set ‘true’ after some point. This is useful
   when we multiply the vector with lower triangular matrix.

   ‘true’ means that a parameter was changed. */

class ParameterSignal
{
protected:
  std::vector<bool> data;
public:
  ParameterSignal(int n);
  ParameterSignal(const ParameterSignal &sig) = default;
  ~ParameterSignal() = default;
  void signalAfter(int l);
  bool
  operator[](int i) const
  {
    return data[i];
  }
  std::vector<bool>::reference
  operator[](int i)
  {
    return data[i];
  }
};

/* This is the abstract class for vector function. At this level of abstraction
   we only need to know size of input vector and a size of output vector.

   The important thing here is a clone method, we will need to make hard copies
   of vector functions since the evaluations are not const. The hardcopies
   apply for parallelization. */

class VectorFunction
{
protected:
  int in_dim;
  int out_dim;
public:
  VectorFunction(int idim, int odim)
    : in_dim(idim), out_dim(odim)
  {
  }
  VectorFunction(const VectorFunction &func) = default;
  virtual ~VectorFunction() = default;
  virtual std::unique_ptr<VectorFunction> clone() const = 0;
  virtual void eval(const Vector &point, const ParameterSignal &sig, Vector &out) = 0;
  int
  indim() const
  {
    return in_dim;
  }
  int
  outdim() const
  {
    return out_dim;
  }
};

/* This makes ‘n’ copies of VectorFunction. The first constructor make exactly
   ‘n’ new copies, the second constructor copies only the pointer to the first
   and others are hard (real) copies.

   The class is useful for making a given number of copies at once, and this
   set can be reused many times if we need mupliple copis of the function (for
   example for paralelizing the code). */

class VectorFunctionSet
{
private:
  // Stores the hard copies made by the class
  std::vector<std::unique_ptr<VectorFunction>> func_copies;
protected:
  std::vector<VectorFunction *> funcs;
public:
  VectorFunctionSet(const VectorFunction &f, int n);
  VectorFunctionSet(VectorFunction &f, int n);
  ~VectorFunctionSet() = default;
  VectorFunction &
  getFunc(int i)
  {
    return *(funcs[i]);
  }
  int
  getNum() const
  {
    return funcs.size();
  }
};

/* This class wraps another VectorFunction to allow integration of a function
   through normally distributed inputs. Namely, if one wants to integrate

         1
    ─────────── ∫ f(x)e^{−½xᵀ|Σ|⁻¹x}dx
    √{(2π)ⁿ|Σ|}

   then if we write Σ=AAᵀ and x=√2·Ay, we get integral

         1                                         1                                     
    ─────────── ∫ f(√2·Ay)e^{−½yᵀy} √(2ⁿ)|A|dy = ───── ∫ f(√2·Ay)e^{−½yᵀy}dy
    √{(2π)ⁿ|Σ|}                                  √(πⁿ)                               

   which means that a given function f we have to wrap to yield a function

    g(y) = 1/√(πⁿ) f(√2·Ay).

   This is exactly what this class is doing. This transformation is useful
   since the Gauss-Hermite points and weights are defined for weighting
   function e^{−y²}, so this transformation allows using Gauss-Hermite
   quadratures seemlessly in a context of integration through normally
   distributed inputs.

   The class maintains a pointer to the function f. When the object is
   constructed by the first constructor, the f is assumed to be owned by the
   caller. If the object of this class is copied, then f is copied and hence
   stored in a std::unique_ptr. The second constructor takes a smart pointer to
   the function and in that case the class takes ownership of f. */

class GaussConverterFunction : public VectorFunction
{
private:
  std::unique_ptr<VectorFunction> func_storage;
protected:
  VectorFunction *func;
  GeneralMatrix A;
  double multiplier;
public:
  GaussConverterFunction(VectorFunction &f, const GeneralMatrix &vcov);
  GaussConverterFunction(std::unique_ptr<VectorFunction> f, const GeneralMatrix &vcov);
  GaussConverterFunction(const GaussConverterFunction &f);
  ~GaussConverterFunction() override = default;
  std::unique_ptr<VectorFunction>
  clone() const override
  {
    return std::make_unique<GaussConverterFunction>(*this);
  }
  void eval(const Vector &point, const ParameterSignal &sig, Vector &out) override;
private:
  double calcMultiplier() const;
  void calcCholeskyFactor(const GeneralMatrix &vcov);
};

#endif
