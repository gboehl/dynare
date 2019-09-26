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
 * along with Dynare.  If not, see <http://www.gnu.org/licenses/>.
 */

#include "quadrature.hh"
#include "precalc_quadrature.hh"

#include <cmath>

void
OneDPrecalcQuadrature::calcOffsets()
{
  offsets[0] = 0;
  for (int i = 1; i < num_levels; i++)
    offsets[i] = offsets[i-1] + num_points[i-1];
}

GaussHermite::GaussHermite()
  : OneDPrecalcQuadrature(gh_num_levels, gh_num_points, gh_weights, gh_points)
{
}

GaussLegendre::GaussLegendre()
  : OneDPrecalcQuadrature(gl_num_levels, gl_num_points, gl_weights, gl_points)
{
}
