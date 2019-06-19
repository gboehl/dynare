/*
 * Copyright © 2006 Ondra Kamenik
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
 * along with Dynare.  If not, see <http://www.gnu.org/licenses/>.
 */

#ifndef OGU_NLSOLVE_H
#define OGU_NLSOLVE_H

#include "twod_matrix.hh"
#include "journal.hh"

#include <cmath>

namespace ogu
{
  class OneDFunction
  {
  public:
    virtual ~OneDFunction() = default;
    virtual double eval(double) = 0;
  };

  class GoldenSectionSearch
  {
  protected:
    constexpr static double tol = 1.e-4;
    /* This is equal to the golden section ratio. */
    static double golden;
  public:
    static double search(OneDFunction &f, double x1, double x2);
  protected:
    /* This initializes a bracket by moving x2 and b (as a golden section of
       x1,x2) so that f(x1)>f(b) && f(b)<f(x2). The point x1 is not moved,
       since it is considered as reliable and f(x1) is supposed to be finite.
       If initialization of a bracket succeeded, then [x1,b,x2] is the bracket
       and true is returned. Otherwise, b is the minimum found and false is
       returned. */
    static bool init_bracket(OneDFunction &f, double x1, double &x2, double &b);
    /* This supposes that f(x1) is finite and it moves x2 toward x1 until x2
       and b (as a golden section of x1,x2) are finite. If succeeded, the
       routine returns true and x2, and b. Otherwise, it returns false. */
    static bool search_for_finite(OneDFunction &f, double x1, double &x2, double &b);
  };

  class VectorFunction
  {
  public:
    VectorFunction() = default;
    virtual ~VectorFunction() = default;
    virtual int inDim() const = 0;
    virtual int outDim() const = 0;
    /* Check dimensions of eval parameters. */
    void check_for_eval(const ConstVector &in, Vector &out) const;
    /* Evaluate the vector function. */
    virtual void eval(const ConstVector &in, Vector &out) = 0;
  };

  class Jacobian : public TwoDMatrix
  {
  public:
    Jacobian(int n) : TwoDMatrix(n, n)
    {
    }
    ~Jacobian() override = default;
    virtual void eval(const Vector &in) = 0;
  };

  class NLSolver : public OneDFunction
  {
  protected:
    Journal &journal;
    VectorFunction &func;
    Jacobian &jacob;
    const int max_iter;
    const double tol;
  private:
    Vector xnewton;
    Vector xcauchy;
    Vector x;
  public:
    NLSolver(VectorFunction &f, Jacobian &j, int maxit, double tl, Journal &jr)
      : journal(jr), func(f), jacob(j), max_iter(maxit), tol(tl),
        xnewton(f.inDim()), xcauchy(f.inDim()), x(f.inDim())
    {
      xnewton.zeros(); xcauchy.zeros(); x.zeros();
    }
    ~NLSolver() override = default;
    /* Returns true if the problem has converged. xx as input is the starting
       value, as output it is a solution. */
    bool solve(Vector &xx, int &iter);
    /* To implement OneDFunction interface. It returns func(xx)ᵀ·func(xx),
       where xx=x+lambda·xcauchy+(1−lambda)·xnewton. It is non-const only
       because it calls func, x, xnewton, xcauchy is not changed. */
    double eval(double lambda) override;
  };
};

#endif
