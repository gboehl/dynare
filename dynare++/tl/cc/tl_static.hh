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
 * along with Dynare.  If not, see <http://www.gnu.org/licenses/>.
 */

// Tensor library static data.

/* The purpose of this file is to make a unique static variable which
   would contain all other static variables and be responsible for their
   correct initialization and destruction. The variables include an
   equivalence bundle and a permutation bundle. Both depend on dimension of the
   problem, and maximum number of variables.

   TLStatic::init() must be called at the beginning of the program, as
   soon as dimension and number of variables is known. */

#ifndef TL_STATIC_H
#define TL_STATIC_H

#include "equivalence.hh"
#include "permutation.hh"

namespace TLStatic
{
  const EquivalenceSet &getEquiv(int n);
  const PermutationSet &getPerm(int n);
  void init(int dim, int nvar);
};

#endif
