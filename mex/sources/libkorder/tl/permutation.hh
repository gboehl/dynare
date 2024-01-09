/*
 * Copyright © 2004 Ondra Kamenik
 * Copyright © 2019-2024 Dynare Team
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

// Permutations.

/* The permutation class is useful when describing a permutation of
   indices in permuted symmetry tensor. This tensor comes to existence,
   for instance, as a result of the following tensor multiplication:

    [g_y³]_γ₁γ₂γ₃ [g_yu]^γ₁_α₁β₃ [g_yu]^γ₂_α₂β₁ [g_u]^γ₃_β₂

   If this operation is done by a Kronecker product of unfolded tensors, the
   resulting tensor has permuted indices. So, in this case the permutation is
   implied by the equivalence: { {0,4}, {1,3}, {2} }. This results in a
   permutation which maps indices (0,1,2,3,4)↦(0,2,4,3,1).

   The other application of Permutation class is to permute indices with the
   same permutation as done during sorting.

   Here we only define an abstraction for the permutation defined by an
   equivalence. Its basic operation is to apply the permutation to the integer
   sequence. The application is right (or inner), in sense that it works on
   indices of the sequence not items of the sequence. More formally s∘m ≠ m∘s.
   In here, the application of the permutation defined by map m is s∘m.

   Also, we need PermutationSet class which contains all permutations
   of an n-element set, and a bundle of permutations PermutationBundle
   which contains all permutation sets up to a given number.
*/

#ifndef PERMUTATION_HH
#define PERMUTATION_HH

#include "equivalence.hh"
#include "int_sequence.hh"

#include <vector>

/* The permutation object will have a map, which defines mapping of indices
   (0,1,…,n-1)↦(m₀,m₁,…, mₙ₋₁). The map is the sequence (m₀,m₁,…, mₙ₋₁). When
   the permutation with the map m is applied on sequence s, it permutes its
   indices: s∘id↦s∘m.

   So we have one constructor from equivalence, then a method apply(),
   and finally a method tailIdentity() which returns a number of trailing
   indices which yield identity. Also we have a constructor calculating
   map, which corresponds to permutation in sort. This is, we want
   (sorted s)∘m = s.
*/

class Permutation
{
protected:
  IntSequence permap;

public:
  explicit Permutation(int len) : permap(len)
  {
    for (int i = 0; i < len; i++)
      permap[i] = i;
  }
  explicit Permutation(const Equivalence& e) : permap(e.getN())
  {
    e.trace(permap);
  }
  Permutation(const Equivalence& e, const Permutation& per) : permap(e.getN())
  {
    e.trace(permap, per);
  }
  explicit Permutation(const IntSequence& s) : permap(s.size())
  {
    computeSortingMap(s);
  };
  Permutation(const Permutation& p1, const Permutation& p2) : permap(p2.permap)
  {
    p1.apply(permap);
  }
  Permutation(const Permutation& p, int i) : permap(p.permap.insert(p.size(), i))
  {
  }
  [[nodiscard]] bool
  operator==(const Permutation& p) const
  {
    return permap == p.permap;
  }
  [[nodiscard]] int
  size() const
  {
    return permap.size();
  }
  void
  print() const
  {
    permap.print();
  }
  void apply(const IntSequence& src, IntSequence& tar) const;
  void apply(IntSequence& tar) const;
  void inverse();
  [[nodiscard]] int tailIdentity() const;
  [[nodiscard]] const IntSequence&
  getMap() const
  {
    return permap;
  }
  IntSequence&
  getMap()
  {
    return permap;
  }

protected:
  void computeSortingMap(const IntSequence& s);
};

/* The PermutationSet maintains an array of of all permutations. The default
   constructor constructs one element permutation set of one element sets. The
   second constructor constructs a new permutation set over n from all
   permutations over n-1. The parameter n needs not to be provided, but it
   serves to distinguish the constructor from copy constructor.

   The method getPreserving() returns a factor subgroup of permutations, which
   are invariants with respect to the given sequence. This are all permutations
   p yielding p∘s = s, where s is the given sequence. */

class PermutationSet
{
  int order {1};
  int size {1};
  std::vector<Permutation> pers;

public:
  PermutationSet();
  PermutationSet(const PermutationSet& ps, int n);
  [[nodiscard]] int
  getNum() const
  {
    return size;
  }
  [[nodiscard]] const Permutation&
  get(int i) const
  {
    return pers[i];
  }
  [[nodiscard]] std::vector<Permutation> getPreserving(const IntSequence& s) const;
};

/* The permutation bundle encapsulates all permutations sets up to some
   given dimension. */

class PermutationBundle
{
  std::vector<PermutationSet> bundle;

public:
  PermutationBundle(int nmax);
  [[nodiscard]] const PermutationSet& get(int n) const;
  void generateUpTo(int nmax);
};

#endif
