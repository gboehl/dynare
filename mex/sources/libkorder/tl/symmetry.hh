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

// Symmetry.

/* Symmetry is an abstraction for a term of the form y³u². It manages only
   indices, not the variable names. So if one uses this abstraction, it must
   be kept in mind that y is the first and u is the second.

   In fact, the symmetry is a special case of equivalence, but its
   implementation is much simpler. We do not need an abstraction for the
   term “yyuyu” but due to Green theorem we can have term y³u². That
   is why the equivalence is too general for our purposes.

   One of a main purposes of the tensor library is to calculate something like:

                                                ₗ
    [B_y²u³]_α₁α₂β₁β₂β₃ = [f_zˡ]_γ₁…γₗ    ∑     ∏  [g_{s^|cₘ|}]_cₘ(α,β)^γₘ
                                       c∈ℳₗ,₅ ᵐ⁼¹

   We must be able to calculate a symmetry induced by symmetry y²u³ and by an
   equivalence class from equivalence c. For equivalence class {0,4} the
   induced symmetry is “yu”, since we pick first and fifth variable from y²u³.
   For a given outer symmetry, the class InducedSymmetries does this for all
   classes of a given equivalence.

   We need also to cycle through all possible symmetries yielding the
   given dimension. For this purpose we define classes SymmetrySet and
   symiterator.

   The symmetry is implemented as IntSequence, in fact, it inherits
   from it. */

#ifndef SYMMETRY_HH
#define SYMMETRY_HH

#include "equivalence.hh"
#include "int_sequence.hh"

#include <initializer_list>
#include <list>
#include <memory>
#include <utility>
#include <vector>

/* Clear. The method isFull() returns true if and only if the symmetry
   allows for any permutation of indices.

   WARNING: Symmetry(n) and Symmetry{n} are not the same. The former
   initializes a symmetry of n elements, while the latter is a full symmetry of
   order n. This is similar to the behaviour of std::vector. */

class Symmetry : public IntSequence
{
public:
  // Constructor allocating a given length of (zero-initialized) data
  explicit Symmetry(int len) : IntSequence(len, 0)
  {
  }
  /* Constructor using an initializer list, that gives the contents of the
     Symmetry. Typically used for symmetries of the form yⁿ, yⁿuᵐ, yⁿuᵐσᵏ */
  Symmetry(std::initializer_list<int> init) : IntSequence(std::move(init))
  {
  }
  // Constructor of implied symmetry for a symmetry and an equivalence class
  Symmetry(const Symmetry& s, const OrdSequence& cl);
  /* Subsymmetry, which takes the given length of symmetry from the end (shares
     data pointer) */
  Symmetry(Symmetry& s, int len) : IntSequence(s, s.size() - len, s.size())
  {
  }

  [[nodiscard]] int
  num() const
  {
    return size();
  }
  [[nodiscard]] int
  dimen() const
  {
    return sum();
  }
  [[nodiscard]] int findClass(int i) const;
  [[nodiscard]] bool isFull() const;
};

/* This is an iterator that iterates over all symmetries of given length and
   dimension (see the SymmetrySet class for details).

   The beginning iterator is (0, …, 0, dim).
   Increasing it gives (0, … , 1, dim−1)
   The just-before-end iterator is (dim, 0, …, 0)
   The past-the-end iterator is (dim+1, 0, …, 0)

   The constructor creates the iterator which starts from the given symmetry
   symmetry (beginning). */

class symiterator
{
  const int dim;
  Symmetry run;

public:
  symiterator(int dim_arg, Symmetry run_arg);
  ~symiterator() = default;
  symiterator& operator++();
  const Symmetry&
  operator*() const
  {
    return run;
  }
  bool
  operator==(const symiterator& it)
  {
    return dim == it.dim && run == it.run;
  }
  bool
  operator!=(const symiterator& it)
  {
    return !operator==(it);
  }
};

/* The class SymmetrySet defines a set of symmetries of the given length
   having given dimension (i.e. it represents all the lists of integers of
   length ‘len’ and of sum equal to ‘dim’). It does not store all the
   symmetries, it is just a convenience class for iteration.

   The typical usage of the abstractions for SymmetrySet and
   symiterator is as follows:

     for (auto &si : SymmetrySet(6, 4))

   It goes through all symmetries of lenght 4 having dimension 6. One can use
   ‘si’ as the symmetry in the body. */

class SymmetrySet
{
public:
  const int len;
  const int dim;
  SymmetrySet(int dim_arg, int len_arg) : len(len_arg), dim(dim_arg)
  {
  }
  [[nodiscard]] symiterator
  begin() const
  {
    Symmetry run(len);
    run[len - 1] = dim;
    return {dim, run};
  }
  [[nodiscard]] symiterator
  end() const
  {
    Symmetry run(len);
    run[0] = dim + 1;
    return {dim, run};
  }
};

/* This simple abstraction just constructs a vector of induced
   symmetries from the given equivalence and outer symmetry. A
   permutation might optionally permute the classes of the equivalence. */

class InducedSymmetries : public std::vector<Symmetry>
{
public:
  InducedSymmetries(const Equivalence& e, const Symmetry& s);
  InducedSymmetries(const Equivalence& e, const Permutation& p, const Symmetry& s);
  void print() const;
};

#endif
