/*
 * Copyright © 2021-2024 Dynare Team
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

#ifndef OBJECTIVE_M_HH
#define OBJECTIVE_M_HH

#include <vector>

#include <dynmex.h>

#include "twod_matrix.hh"

// Handles calls to <model>/+objective/static.m
class ObjectiveMFile
{
private:
  int ntt; // Size of vector of temporary terms
  const std::string ObjectiveMFilename;
  static void unpackSparseMatrixAndCopyIntoTwoDMatData(mxArray* sparseMat, TwoDMatrix& tdm);

public:
  explicit ObjectiveMFile(const std::string& modName, int ntt_arg);
  void eval(const Vector& y, const Vector& x, const Vector& params, Vector& residual,
            std::vector<TwoDMatrix>& md);
};
#endif
