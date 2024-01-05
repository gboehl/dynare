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

// Integer sequence.

/* Here we define an auxiliary abstraction for a sequence of integers. The
   basic functionality is to hold an ordered sequence of integers with
   constant length. We prefer using this simple class rather than
   std::vector<int> since it is more efficient for our purposes.

   The class is used in index of a tensor, in symmetry definition, in
   Kronecker product dimensions, or as a class of an equivalence. The
   latter case is not ordered, but we always order equivalence classes in
   order to ensure unique representativeness. For almost all cases we
   need the integer sequence to be ordered (sort), or monotonize (indices
   of folded tensors), or partially monotonize (indices of folded tensors
   not fully symmetric), or calculate a product of all members or only of
   a part (used in Kronecker product dimensions). When we calculate
   offsets in folded tensors, we need to obtain a number of the same
   items in the front (getPrefixLength()), and also to add some integer
   number to all items.

   Also, we need to construct a subsequence of a sequence. */

#ifndef INT_SEQUENCE_HH
#define INT_SEQUENCE_HH

#include <algorithm>
#include <initializer_list>
#include <utility>
#include <vector>

/* The implementation of IntSequence is straightforward. It has a pointer
   ‘data’, an ‘offset’ integer indicating the beginning of the data relatively
   to the pointer and a ‘length’ of the sequence.

   WARNING: IntSequence(n) and IntSequence{n} are not the same (parentheses
   versus braces). The former initializes a sequence of length n, while the
   latter constructs a sequence of a single element equal to n. This is similar
   to the behaviour of std::vector. */

class Symmetry;
class IntSequence
{
  int* data;
  int length;
  bool destroy {true};

public:
  // Constructor allocating a given length of (uninitialized) data
  explicit IntSequence(int l) : data {new int[l]}, length {l}
  {
  }
  // Constructor allocating and then initializing all members to a given number
  IntSequence(int l, int n) : data {new int[l]}, length {l}
  {
    std::fill_n(data, length, n);
  }
  /* Constructor using an initializer list (gives the contents of the
     IntSequence, similarly to std::vector) */
  IntSequence(std::initializer_list<int> init) :
      data {new int[init.size()]}, length {static_cast<int>(init.size())}
  {
    std::copy(init.begin(), init.end(), data);
  }
  // Copy constructor
  IntSequence(const IntSequence& s) : data {new int[s.length]}, length {s.length}
  {
    std::copy_n(s.data, length, data);
  }
  // Move constructor
  IntSequence(IntSequence&& s) noexcept :
      data {std::exchange(s.data, nullptr)},
      length {std::exchange(s.length, 0)},
      destroy {std::exchange(s.destroy, false)}
  {
  }
  // Subsequence constructor (which shares the data pointer)
  IntSequence(IntSequence& s, int i1, int i2) :
      data {s.data + i1}, length {i2 - i1}, destroy {false}
  {
  }
  // Subsequence constructor (without pointer sharing)
  IntSequence(const IntSequence& s, int i1, int i2) : data {new int[i2 - i1]}, length {i2 - i1}
  {
    std::copy_n(s.data + i1, length, data);
  }
  /* Unfolds a given integer sequence with respect to a given symmetry. If for
     example the sequence is (a,b) and the symmetry is (2,3), then the
     result is (a,a,b,b,b). */
  [[nodiscard]] IntSequence unfold(const Symmetry& sy) const;

  /* Constructs a symmetry from the integer sequence (supposed to be ordered) as
     a symmetry counting successively equal items. For instance the sequence
     (a,a,a,b,c,c,d,d,d,d) produces symmetry (3,1,2,4). */
  [[nodiscard]] Symmetry getSymmetry() const;

  IntSequence& operator=(const IntSequence& s);
  IntSequence& operator=(IntSequence&& s);
  virtual ~IntSequence()
  {
    if (destroy)
      delete[] data;
  }
  bool operator==(const IntSequence& s) const;
  int&
  operator[](int i)
  {
    return data[i];
  }
  int
  operator[](int i) const
  {
    return data[i];
  }
  [[nodiscard]] int
  size() const
  {
    return length;
  }

  /* We provide two orderings. The first operator<() is the linear
     lexicographic ordering, the second less() is the non-linear Cartesian
     ordering. */
  bool operator<(const IntSequence& s) const;
  bool
  operator<=(const IntSequence& s) const
  {
    return (operator==(s) || operator<(s));
  }
  [[nodiscard]] bool lessEq(const IntSequence& s) const;
  [[nodiscard]] bool less(const IntSequence& s) const;

  // Inserts an element into an ordered sequence
  [[nodiscard]] IntSequence insert(int i) const;
  // Inserts an element at a given position
  /* For appending at the end, use pos = size() */
  [[nodiscard]] IntSequence insert(int i, int pos) const;

  // In-place sort the sequence in increasing order
  void sort();

  /* Monotonize the sequence: if an item is less then its predecessor, it is
     equalized. */
  void monotone();

  /* Partially monotonize the sequence. The partitioning is done by a
   symmetry. So the subsequence given by the symmetry classes are
   monotonized. For example, if the symmetry is y²u³, and the
   IntSequence is (5,3,1,6,4), the result is (5,5,1,6,6). */
  void pmonotone(const Symmetry& s);

  // Returns the sum of all elements. Useful for symmetries
  [[nodiscard]] int sum() const;

  /* This returns product of elements between indices i1 (included) and i2
   (excluded). Useful for Kronecker product dimensions. */
  [[nodiscard]] int mult(int i1, int i2) const;

  // Returns the product of all elements
  [[nodiscard]] int
  mult() const
  {
    return mult(0, length);
  }
  void add(int i);
  void add(int f, const IntSequence& s);

  /* Return the number of identical elements at the beginning of the sequence. */
  [[nodiscard]] int getPrefixLength() const;

  /* This returns a number of distinct items in the sequence. It supposes that
     the sequence is ordered. Returns zero on the empty sequence. */
  [[nodiscard]] int getNumDistinct() const;

  /* This returns a maximum of the sequence. If the sequence is empty, it
     returns the least possible int value. */
  [[nodiscard]] int getMax() const;

  [[nodiscard]] bool isPositive() const;
  [[nodiscard]] bool isConstant() const;
  [[nodiscard]] bool isSorted() const;

  void print() const;
  /*                                   ⎛sum(this)⎞
     Computes multinomial coefficient: ⎝  this   ⎠
     (where the lower line consists of the sequence of integers stored by ‘this’)

     WARNING: this operation is destructive; make a copy if you want to keep
     the original sequence */
  int noverseq();
};

#endif
