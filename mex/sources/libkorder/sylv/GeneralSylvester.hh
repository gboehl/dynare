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

#ifndef GENERAL_SYLVESTER_H
#define GENERAL_SYLVESTER_H

#include "SimilarityDecomp.hh"
#include "SylvMatrix.hh"
#include "SylvesterSolver.hh"

#include <memory>

class GeneralSylvester
{
  SylvParams pars;
  int order;
  const SqSylvMatrix a;
  const SylvMatrix b;
  const SqSylvMatrix c;
  SylvMatrix d;
  bool solved;
  std::unique_ptr<SchurDecompZero> bdecomp;
  std::unique_ptr<SimilarityDecomp> cdecomp;
  std::unique_ptr<SylvesterSolver> sylv;

public:
  // Construct with my copy of d
  GeneralSylvester(int ord, int n, int m, int zero_cols, const ConstVector& da,
                   const ConstVector& db, const ConstVector& dc, const ConstVector& dd,
                   const SylvParams& ps);
  GeneralSylvester(int ord, int n, int m, int zero_cols, const ConstVector& da,
                   const ConstVector& db, const ConstVector& dc, const ConstVector& dd,
                   bool alloc_for_check = false);
  // Construct with provided storage for d
  GeneralSylvester(int ord, int n, int m, int zero_cols, const ConstVector& da,
                   const ConstVector& db, const ConstVector& dc, Vector& dd,
                   bool alloc_for_check = false);
  GeneralSylvester(int ord, int n, int m, int zero_cols, const ConstVector& da,
                   const ConstVector& db, const ConstVector& dc, Vector& dd, const SylvParams& ps);
  virtual ~GeneralSylvester() = default;
  int
  getM() const
  {
    return c.nrows();
  }
  int
  getN() const
  {
    return a.nrows();
  }
  const double*
  getResult() const
  {
    return d.base();
  }
  const SylvParams&
  getParams() const
  {
    return pars;
  }
  SylvParams&
  getParams()
  {
    return pars;
  }
  void solve();
  void check(const ConstVector& ds);

private:
  void init();
};

#endif /* GENERAL_SYLVESTER_H */
