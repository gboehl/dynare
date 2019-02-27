// Copyright 2004, Ondra Kamenik

#include "int_sequence.hh"
#include "symmetry.hh"
#include "tl_exception.hh"
#include "pascal_triangle.hh"

#include <iostream>
#include <limits>
#include <numeric>

/* This unfolds a given integer sequence with respect to the given
   symmetry. If for example the symmetry is $(2,3)$, and the sequence is
   $(a,b)$, then the result is $(a,a,b,b,b)$. */

IntSequence::IntSequence(const Symmetry &sy, const IntSequence &se)
  : data{new int[sy.dimen()]}, length{sy.dimen()}
{
  int k = 0;
  for (int i = 0; i < sy.num(); i++)
    for (int j = 0; j < sy[i]; j++, k++)
      operator[](k) = se[i];
}

/* This constructs an implied symmetry (implemented as |IntSequence|
   from a more general symmetry and equivalence class (implemented as
   |vector<int>|). For example, let the general symmetry be $y^3u^2$ and
   the equivalence class is $\{0,4\}$ picking up first and fifth
   variable, we calculate symmetry (at this point only |IntSequence|)
   corresponding to the picked variables. These are $yu$. Thus the
   constructed sequence must be $(1,1)$, meaning that we picked one $y$
   and one $u$. */

IntSequence::IntSequence(const Symmetry &sy, const std::vector<int> &se)
  : data{new int[sy.num()]}, length{sy.num()}
{
  TL_RAISE_IF(sy.dimen() <= se[se.size()-1],
              "Sequence is not reachable by symmetry in IntSequence()");
  for (int i = 0; i < length; i++)
    operator[](i) = 0;

  for (int i : se)
    operator[](sy.findClass(i))++;
}

/* This constructs an ordered integer sequence from the given ordered
   sequence inserting the given number to the sequence. */

IntSequence
IntSequence::insert(int i) const
{
  IntSequence r(size()+1);
  int j;
  for (j = 0; j < size() && operator[](j) < i; j++)
    r[j] = operator[](j);
  r[j] = i;
  for (; j < size(); j++)
    r[j+1] = operator[](j);
  return r;
}

IntSequence
IntSequence::insert(int i, int pos) const
{
  TL_RAISE_IF(pos < 0 || pos > size(),
              "Wrong position for IntSequence::insert()");
  IntSequence r(size()+1);
  int j;
  for (j = 0; j < pos; j++)
    r[j] = operator[](j);
  r[j] = i;
  for (; j < size(); j++)
    r[j+1] = operator[](j);
  return r;
}

IntSequence &
IntSequence::operator=(const IntSequence &s)
{
  TL_RAISE_IF(length != s.length, "Wrong length for in-place IntSequence::operator=");
  std::copy_n(s.data, length, data);
  return *this;
}

IntSequence &
IntSequence::operator=(IntSequence &&s)
{
  TL_RAISE_IF(length != s.length, "Wrong length for in-place IntSequence::operator=");
  std::copy_n(s.data, length, data);
  return *this;
}

bool
IntSequence::operator==(const IntSequence &s) const
{
  return std::equal(data, data+length,
                    s.data, s.data+s.length);
}

bool
IntSequence::operator<(const IntSequence &s) const
{
  return std::lexicographical_compare(data, data+length,
                                      s.data, s.data+s.length);
}

bool
IntSequence::lessEq(const IntSequence &s) const
{
  TL_RAISE_IF(size() != s.size(),
              "Sequence with different lengths in IntSequence::lessEq");

  int i = 0;
  while (i < size() && operator[](i) <= s[i])
    i++;
  return (i == size());
}

bool
IntSequence::less(const IntSequence &s) const
{
  TL_RAISE_IF(size() != s.size(),
              "Sequence with different lengths in IntSequence::less");

  int i = 0;
  while (i < size() && operator[](i) < s[i])
    i++;
  return (i == size());
}

void
IntSequence::sort()
{
  std::sort(data, data+length);
}

/* Here we monotonize the sequence. If an item is less then its
   predecessor, it is equalized. */

void
IntSequence::monotone()
{
  for (int i = 1; i < length; i++)
    if (operator[](i-1) > operator[](i))
      operator[](i) = operator[](i-1);
}

/* This partially monotones the sequence. The partitioning is done by a
   symmetry. So the subsequence given by the symmetry classes are
   monotonized. For example, if the symmetry is $y^2u^3$, and the
   |IntSequence| is $(5,3,1,6,4)$, the result is $(5,5,1,6,6)$. */

void
IntSequence::pmonotone(const Symmetry &s)
{
  int cum = 0;
  for (int i = 0; i < s.num(); i++)
    {
      for (int j = cum + 1; j < cum + s[i]; j++)
        if (operator[](j-1) > operator[](j))
          operator[](j) = operator[](j-1);
      cum += s[i];
    }
}

/* This returns sum of all elements. Useful for symmetries. */

int
IntSequence::sum() const
{
  return std::accumulate(data, data+length, 0);
}

/* This returns product of subsequent items. Useful for Kronecker product
   dimensions. */

int
IntSequence::mult(int i1, int i2) const
{
  return std::accumulate(data+i1, data+i2,
                         1, std::multiplies<int>());
}

/* Return a number of the same items in the beginning of the sequence. */

int
IntSequence::getPrefixLength() const
{
  int i = 0;
  while (i+1 < size() && operator[](i+1) == operator[](0))
    i++;
  return i+1;
}

/* This returns a number of distinct items in the sequence. It supposes
   that the sequence is ordered. For the empty sequence it returns zero. */

int
IntSequence::getNumDistinct() const
{
  int res = 0;
  if (length > 0)
    res++;
  for (int i = 1; i < length; i++)
    if (operator[](i) != operator[](i-1))
      res++;
  return res;
}

/* This returns a maximum of the sequence. If the sequence is empty, it
   returns the least possible |int| value. */

int
IntSequence::getMax() const
{
  if (length == 0)
    return std::numeric_limits<int>::min();
  return *std::max_element(data, data+length);
}

void
IntSequence::add(int i)
{
  for (int j = 0; j < size(); j++)
    operator[](j) += i;
}

void
IntSequence::add(int f, const IntSequence &s)
{
  TL_RAISE_IF(size() != s.size(),
              "Wrong sequence length in IntSequence::add");
  for (int j = 0; j < size(); j++)
    operator[](j) += f*s[j];
}

bool
IntSequence::isPositive() const
{
  return std::all_of(data, data+length,
                     [](int x) { return x >= 0; });
}

bool
IntSequence::isConstant() const
{
  if (length < 2)
    return true;
  return std::all_of(data+1, data+length,
                     [this](int x) { return x == operator[](0); });
}

bool
IntSequence::isSorted() const
{
  return std::is_sorted(data, data+length);
}

/* Debug print. */

void
IntSequence::print() const
{
  std::cout << '[';
  for (int i = 0; i < size(); i++)
    std::cout << operator[](i) << ' ';
  std::cout << ']' << std::endl;
}

/* Here we calculate the multinomial coefficients
    $\left(\matrix{a\cr b_1,\ldots,b_n}\right)$, where $a=b_1+\ldots+b_n$.

   See:
    https://en.wikipedia.org/wiki/Binomial_coefficient#Generalization_to_multinomials
    https://en.wikipedia.org/wiki/Multinomial_theorem

   For n=1, the coefficient is equal to 1.
   For n=2, the multinomial coeffs correspond to the binomial coeffs, i.e. the binomial
    (a; b) is equal to the multinomial (a; b,a-b).
   For n>=3, we have the identity
   $$\left(\matrix{a\cr b_1,\ldots,b_n}\right)=\left(\matrix{b_1+b_2\cr b_1}\right)\cdot
   \left(\matrix{a\cr b_1+b_2,b_3,\ldots,b_n}\right)$$ (where the first factor
   on the right hand side is to be interpreted as a binomial coefficient)

   This number is exactly a number of unfolded indices corresponding to one
   folded index, where the sequence $b_1,\ldots,b_n$ is the symmetry of the
   index. This can be easily seen if the multinomial coefficient is interpreted
   as the number of unique permutations of a word, where $a$ is the length of
   the word, $n$ is the number of distinct letters, and the $b_i$ are the
   number of repetitions of each letter. For example, for a symmetry of the
   form $y^4 u^2 v^3$, we want to compute the number of permutations of the word
   $yyyyuuvvv$. This is equal to the multinomial coefficient (9; 4,2,3). */

int
IntSequence::noverseq()
{
  if (size() == 0 || size() == 1)
    return 1;
  data[1] += data[0];
  return PascalTriangle::noverk(data[1], data[0]) * IntSequence(*this, 1, size()).noverseq();
}
