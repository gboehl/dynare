/*
 * Copyright (C) 2010-2019 Dynare Team
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

#ifndef _DYNAMICMODELAC_HH
#define _DYNAMICMODELAC_HH

#include "twod_matrix.hh"

class DynamicModelAC
{
public:
  virtual ~DynamicModelAC() = default;
  virtual void eval(const Vector &y, const Vector &x, const Vector &params, const Vector &ySteady,
                    Vector &residual, TwoDMatrix *g1, TwoDMatrix *g2, TwoDMatrix *g3) = 0;
};
#endif
