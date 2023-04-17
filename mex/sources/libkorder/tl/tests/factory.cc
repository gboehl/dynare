/*
 * Copyright © 2004 Ondra Kamenik
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

#include "factory.hh"

void
Factory::init(const Symmetry &s, const IntSequence &nvs)
{
  IntSequence sym(s);
  decltype(mtgen)::result_type seed = sym[0];
  seed = 256*seed + nvs[0];
  if (sym.size() > 1)
    seed = 256*seed + sym[1];
  if (nvs.size() > 1)
    seed = 256*seed + nvs[0];
  mtgen.seed(seed);
}

void
Factory::init(int dim, int nv)
{
  decltype(mtgen)::result_type seed = dim;
  seed = 256*seed + nv;
  mtgen.seed(seed);
}

double
Factory::get()
{
  return dis(mtgen)-0.5;
}

void
Factory::fillMatrix(TwoDMatrix &m)
{
  Vector &d = m.getData();
  for (int i = 0; i < d.length(); i++)
    d[i] = get();
}

Vector
Factory::makeVector(int n)
{
  init(n, n*n);

  Vector v(n);
  for (int i = 0; i < n; i++)
    v[i] = get();

  return v;
}
