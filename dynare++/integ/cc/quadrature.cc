// Copyright 2005, Ondra Kamenik

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
