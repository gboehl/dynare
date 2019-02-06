// Copyright 2004, Ondra Kamenik

// Integer sequence.

/* Here we define an auxiliary abstraction for a sequence of integers. The
   basic functionality is to hold an ordered sequence of integers with
   constant length. We prefer using this simple class before STL
   |vector<int>| since it is more efficient for our purposes.

   The class is used in index of a tensor, in symmetry definition, in
   Kronecker product dimensions, or as a class of an equivalence. The
   latter case is not ordered, but we always order equivalence classes in
   order to ensure unique representativeness. For almost all cases we
   need the integer sequence to be ordered (sort), or monotonize (indices
   of folded tensors), or partially monotonize (indices of folded tensors
   not fully symmetric), or calculate a product of all members or only of
   a part (used in Kronecker product dimensions). When we calculate
   offsets in folded tensors, we need to obtain a number of the same
   items in the front (|getPrefixLength|), and also to add some integer
   number to all items.

   Also, we need to construct a subsequence of a sequence. */

#ifndef INT_SEQUENCE_H
#define INT_SEQUENCE_H

#include <vector>
#include <memory>
#include <algorithm>

/* The implementation of |IntSequence| is straightforward. It has a pointer
   |data|, an |offset| integer indicating the beginning of the data relatively
   to the pointer and a |length| of the sequence. */

class Symmetry;
class IntSequence
{
  std::shared_ptr<int> data;
  int length;
  int offset{0};
public:
  /* We have a constructor allocating a given length of data, constructor
     allocating and then initializing all members to a given number, a copy
     constructor, a conversion from |vector<int>|, a subsequence
     constructor, a constructor used for calculating implied symmetry from
     a more general symmetry and one equivalence class (see |Symmetry|
     class). Finally we have a constructor which unfolds a sequence with
     respect to a given symmetry and constructor which inserts a given
     number to the ordered sequence or given number to a given position. */

  explicit IntSequence(int l)
    : data{new int[l], [](int *arr) { delete[] arr; }}, length{l}
  {
  }
  IntSequence(int l, int n)
    : data{new int[l], [](int *arr) { delete[] arr; }}, length{l}
  {
    std::fill_n(data.get(), length, n);
  }
  IntSequence(const IntSequence &s)
    : data{new int[s.length], [](int *arr) { delete[] arr; }}, length{s.length}
  {
    std::copy_n(s.data.get()+s.offset, length, data.get());
  }
  IntSequence(IntSequence &&s) = default;
  IntSequence(IntSequence &s, int i1, int i2)
    : data{s.data}, length{i2-i1}, offset{s.offset+i1}
  {
  }
  IntSequence(const IntSequence &s, int i1, int i2)
    : data{new int[i2-i1], [](int *arr) { delete[] arr; }}, length{i2-i1}
  {
    std::copy_n(s.data.get()+s.offset+i1, length, data.get());
  }
  IntSequence(const Symmetry &sy, const std::vector<int> &se);
  IntSequence(const Symmetry &sy, const IntSequence &se);
  IntSequence(int i, const IntSequence &s);
  IntSequence(int i, const IntSequence &s, int pos);

  const IntSequence &operator=(const IntSequence &s);
  const IntSequence &operator=(IntSequence &&s);
  virtual ~IntSequence() = default;
  bool operator==(const IntSequence &s) const;
  bool
  operator!=(const IntSequence &s) const
  {
    return !operator==(s);
  }
  int &
  operator[](int i)
  {
    return data.get()[offset+i];
  }
  int
  operator[](int i) const
  {
    return data.get()[offset+i];
  }
  int
  size() const
  {
    return length;
  }

  /* We provide two orderings. The first |operator<| is the linear
     lexicographic ordering, the second |less| is the non-linear Cartesian
     ordering. */
  bool operator<(const IntSequence &s) const;
  bool
  operator<=(const IntSequence &s) const
  {
    return (operator==(s) || operator<(s));
  }
  bool lessEq(const IntSequence &s) const;
  bool less(const IntSequence &s) const;

  void sort();
  void monotone();
  void pmonotone(const Symmetry &s);
  int sum() const;
  int mult(int i1, int i2) const;
  int
  mult() const
  {
    return mult(0, length);
  }
  void add(int i);
  void add(int f, const IntSequence &s);
  int getPrefixLength() const;
  int getNumDistinct() const;
  int getMax() const;
  bool isPositive() const;
  bool isConstant() const;
  bool isSorted() const;
  void print() const;
};

#endif
