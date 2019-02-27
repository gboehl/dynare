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
#include <algorithm>
#include <initializer_list>

/* The implementation of |IntSequence| is straightforward. It has a pointer
   |data|, an |offset| integer indicating the beginning of the data relatively
   to the pointer and a |length| of the sequence.

   WARNING: IntSequence(n) and IntSequence{n} are not the same. The former
   initializes a sequence of length n, while the latter constructs a sequence
   of a single element equal to n. This is similar to the behaviour of
   std::vector. */

class Symmetry;
class IntSequence
{
  int *data;
  int length;
  bool destroy{true};
public:
  // Constructor allocating a given length of (uninitialized) data
  explicit IntSequence(int l)
    : data{new int[l]}, length{l}
  {
  }
  // Constructor allocating and then initializing all members to a given number
  IntSequence(int l, int n)
    : data{new int[l]}, length{l}
  {
    std::fill_n(data, length, n);
  }
  /* Constructor using an initializer list (gives the contents of the
     IntSequence, similarly to std::vector) */
  IntSequence(std::initializer_list<int> init)
    : data{new int[init.size()]},
      length{static_cast<int>(init.size())}
  {
    std::copy(init.begin(), init.end(), data);
  }
  // Copy constructor
  IntSequence(const IntSequence &s)
    : data{new int[s.length]}, length{s.length}
  {
    std::copy_n(s.data, length, data);
  }
  // Move constructor
  IntSequence(IntSequence &&s)
    : data{s.data}, length{s.length}, destroy{s.destroy}
  {
    s.data = nullptr;
    s.length = 0;
    s.destroy = false;
  }
  // Subsequence constructor (which shares the data pointer)
  IntSequence(IntSequence &s, int i1, int i2)
    : data{s.data+i1}, length{i2-i1}, destroy{false}
  {
  }
  // Subsequence constructor (without pointer sharing)
  IntSequence(const IntSequence &s, int i1, int i2)
    : data{new int[i2-i1]}, length{i2-i1}
  {
    std::copy_n(s.data+i1, length, data);
  }
  /* Unfolds a given integer sequence with respect to a given symmetry. If for
     example the sequence is $(a,b)$ and the symmetry is $(2,3)$, then the
     result is $(a,a,b,b,b)$. */
  IntSequence unfold(const Symmetry &sy) const;

  IntSequence &operator=(const IntSequence &s);
  IntSequence &operator=(IntSequence &&s);
  virtual ~IntSequence()
  {
    if (destroy)
      delete[] data;
  }
  bool operator==(const IntSequence &s) const;
  bool
  operator!=(const IntSequence &s) const
  {
    return !operator==(s);
  }
  int &
  operator[](int i)
  {
    return data[i];
  }
  int
  operator[](int i) const
  {
    return data[i];
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

  // Inserts an element into an ordered sequence
  IntSequence insert(int i) const;
  // Inserts an element at a given position
  /* For appending at the end, use pos = size() */
  IntSequence insert(int i, int pos) const;

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
  // Compute multinomial coefficient (sum(this); this)
  /* Warning: this operation is destructive; make a copy if you want to keep
     the original sequence */
  int noverseq();
};

#endif
