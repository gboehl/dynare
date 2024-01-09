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

// Equivalences.

/* Here we define an equivalence of a set of integers {0, 1, …, k-1}.
   The purpose is clear, in the tensor library we often iterate
   through all equivalences and sum matrices. We need an abstraction for
   an equivalence class, equivalence and a set of all equivalences.

   The equivalence class (which is basically a set of integers) is here
   implemented as ordered integer sequence. The ordered sequence is not
   implemented via IntSequence, but via vector<int> since we need
   insertions. The equivalence is implemented as an ordered list of
   equivalence classes, and equivalence set is a list of equivalences.

   The ordering of the equivalence classes within an equivalence is very
   important. For instance, if we iterate through equivalences for k=5
   and pickup some equivalence class, say { {0,4}, {1,2}, {3} }, we
   then evaluate something like:

   [B_y²u³]_α₁α₂β₁β₂β₃ = … + [f_z³]_γ₁γ₂γ₃ [z_yu]^γ₁_α₁β₃ [z_yu]^γ₂_α₂β₁ [zᵤ]^γ₃_β₂ + …

   If the tensors are unfolded, we can evaluate this expression as
   f_z³·(z_yu ⊗ z_yu ⊗ zᵤ)·P, where P is a suitable permutation of columns of
   the expressions, which permutes them so that the index
   (α₁,β₃,α₂,β₁,β₂) would go to (α₁,α₂,β₁,β₂,β₃).

   The permutation P can be very ineffective (copying great amount of
   small chunks of data) if the equivalence class ordering is chosen
   badly. However, we do not provide any heuristic minimizing a total
   time spent in all permutations. We choose an ordering which orders the
   classes according to their averages, and according to the smallest
   equivalence class element if the averages are the same. */

#ifndef EQUIVALENCE_HH
#define EQUIVALENCE_HH

#include "int_sequence.hh"

#include <compare>
#include <list>
#include <string>
#include <vector>

/* Here is the abstraction for an equivalence class. We implement it as
   vector<int>. We have a constructor for empty class, copy
   constructor. What is important here is the ordering operator
   operator<=>() and methods for addition of an integer, and addition of
   another sequence. Also we provide method has() which returns true if a
   given integer is contained. */

class OrdSequence
{
  std::vector<int> data;

public:
  OrdSequence() : data()
  {
  }
  [[nodiscard]] bool operator==(const OrdSequence& s) const;
  int operator[](int i) const;
  [[nodiscard]] std::partial_ordering operator<=>(const OrdSequence& s) const;
  [[nodiscard]] const std::vector<int>&
  getData() const
  {
    return data;
  }
  [[nodiscard]] int
  length() const
  {
    return data.size();
  }
  void add(int i);
  void add(const OrdSequence& s);
  [[nodiscard]] bool has(int i) const;
  void print(const std::string& prefix) const;

private:
  [[nodiscard]] double average() const;
};

/* Here is the abstraction for the equivalence. It is a list of
   equivalence classes. Also we remember n, which is a size of
   underlying set {0, 1, …, n-1}.

   Method trace() “prints” the equivalence into the integer sequence. */

class Permutation;
class Equivalence
{
private:
  int n;
  std::list<OrdSequence> classes;

public:
  using const_seqit = std::list<OrdSequence>::const_iterator;
  using seqit = std::list<OrdSequence>::iterator;

  // Constructs { {0}, {1}, …, {n-1} }
  explicit Equivalence(int num);
  // Constructs { {0,1,…,n-1 } }
  Equivalence(int num, const std::string& dummy);
  // Copy constructor plus gluing i1 and i2 in one class
  Equivalence(const Equivalence& e, int i1, int i2);

  [[nodiscard]] bool operator==(const Equivalence& e) const;
  [[nodiscard]] int
  getN() const
  {
    return n;
  }
  [[nodiscard]] int
  numClasses() const
  {
    return classes.size();
  }
  void trace(IntSequence& out, int n) const;
  void
  trace(IntSequence& out) const
  {
    trace(out, numClasses());
  }
  void trace(IntSequence& out, const Permutation& per) const;
  void print(const std::string& prefix) const;
  seqit
  begin()
  {
    return classes.begin();
  }
  [[nodiscard]] const_seqit
  begin() const
  {
    return classes.begin();
  }
  seqit
  end()
  {
    return classes.end();
  }
  [[nodiscard]] const_seqit
  end() const
  {
    return classes.end();
  }
  [[nodiscard]] const_seqit find(int i) const;
  seqit find(int i);

protected:
  /* Here we have find methods. We can find an equivalence class having a
     given number or we can find an equivalence class of a given index within
     the ordering.

     We have also an insert() method which inserts a given class
     according to the class ordering. */
  [[nodiscard]] const_seqit findHaving(int i) const;
  seqit findHaving(int i);
  void insert(const OrdSequence& s);
};

/* The EquivalenceSet is a list of equivalences. The unique
   constructor constructs a set of all equivalences over an n-elements
   set. The equivalences are sorted in the list so that equivalences with
   fewer number of classes are in the end.

   The two methods has() and addParents() are useful in the constructor. */

class EquivalenceSet
{
  int n;
  std::list<Equivalence> equis;

public:
  explicit EquivalenceSet(int num);
  void print(const std::string& prefix) const;
  [[nodiscard]] auto
  begin() const
  {
    return equis.begin();
  }
  [[nodiscard]] auto
  end() const
  {
    return equis.end();
  }

private:
  [[nodiscard]] bool has(const Equivalence& e) const;
  void addParents(const Equivalence& e, std::list<Equivalence>& added);
};

/* The equivalence bundle class only encapsulates EquivalenceSet·s
   from 1 up to a given number. It is able to retrieve the equivalence set
   over n-element set for a given n, and also it can generate some more
   sets on request.

   It is fully responsible for storage needed for EquivalenceSet·s. */

class EquivalenceBundle
{
  std::vector<EquivalenceSet> bundle;

public:
  explicit EquivalenceBundle(int nmax);
  [[nodiscard]] const EquivalenceSet& get(int n) const;
  void generateUpTo(int nmax);
};

#endif
