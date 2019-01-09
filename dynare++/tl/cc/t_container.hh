// Copyright 2004, Ondra Kamenik

// Tensor containers.

/* One of primary purposes of the tensor library is to perform one step
   of the Faa Di Bruno formula:
   $$\left[B_{s^k}\right]_{\alpha_1\ldots\alpha_k}=
   [h_{y^l}]_{\gamma_1\ldots\gamma_l}\sum_{c\in M_{l,k}}
   \prod_{m=1}^l\left[g_{s^{\vert c_m\vert}}\right]^{\gamma_m}_{c_m(\alpha)}
   $$
   where $h_{y^l}$ and $g_{s^i}$ are tensors, $M_{l,k}$ is a set of all
   equivalences with $l$ classes of $k$ element set, $c_m$ is $m$-the
   class of equivalence $c$, and $\vert c_m\vert$ is its
   cardinality. Further, $c_m(\alpha)$ is a sequence of $\alpha$s picked
   by equivalence class $c_m$.

   In order to accomplish this operation, we basically need some storage
   of all tensors of the form $\left[g_{s^i}\right]$. Note that $s$ can
   be compound, for instance $s=[y,u]$. Then we need storage for
   $\left[g_{y^3}\right]$, $\left[g_{y^2u}\right]$,
   $\left[g_{yu^5}\right]$, etc.

   We need an object holding all tensors of the same type. Here type
   means an information, that coordinates of the tensors can be of type
   $y$, or $u$. We will group only tensors, whose symmetry is described
   by |Symmetry| class. These are only $y^2u^3$, not $yuyu^2$. So, we are
   going to define a class which will hold tensors whose symmetries are
   of type |Symmetry| and have the same symmetry length (number of
   different coordinate types). Also, for each symmetry there will be at
   most one tensor.

   The class has two purposes: The first is to provide storage (insert
   and retrieve). The second is to perform the above step of Faa Di Bruno. This is
   going through all equivalences with $l$ classes, perform the tensor
   product and add to the result.

   We define a template class |TensorContainer|. From different
   instantiations of the template class we will inherit to create concrete
   classes, for example container of unfolded general symmetric
   tensors. The one step of the Faa Di Bruno (we call it |multAndAdd|) is
   implemented in the concrete subclasses, because the implementation
   depends on storage. Note even, that |multAndAdd| has not a template
   common declaration. This is because sparse tensor $h$ is multiplied by
   folded tensors $g$ yielding folded tensor $B$, but unfolded tensor $h$
   is multiplied by unfolded tensors $g$ yielding unfolded tensor $B$. */

#ifndef T_CONTAINER_H
#define T_CONTAINER_H

#include "symmetry.hh"
#include "gs_tensor.hh"
#include "tl_exception.hh"
#include "tl_static.hh"
#include "sparse_tensor.hh"
#include "equivalence.hh"
#include "rfs_tensor.hh"
#include "Vector.hh"

#include <map>
#include <string>
#include <sstream>

#include <matio.h>

// |ltsym| predicate
/* We need a predicate on strict weak ordering of
   symmetries. */
struct ltsym
{
  bool
  operator()(const Symmetry &s1, const Symmetry &s2) const
  {
    return s1 < s2;
  }
};

/* Here we define the template class for tensor container. We implement
   it as |stl::map|. It is a unique container, no two tensors with same
   symmetries can coexist. Keys of the map are symmetries, values are
   pointers to tensor. The class is responsible for deallocating all
   tensors. Creation of the tensors is done outside.

   The class has integer |n| as its member. It is a number of different
   coordinate types of all contained tensors. Besides intuitive insert
   and retrieve interface, we define a method |fetchTensors|, which for a
   given symmetry and given equivalence calculates symmetries implied by
   the symmetry and all equivalence classes, and fetches corresponding
   tensors in a vector.

   Also, each instance of the container has a reference to
   |EquivalenceBundle| which allows an access to equivalences. */

template<class _Ttype>
class TensorContainer
{
protected:
  typedef const _Ttype *_const_ptr;
  typedef _Ttype *_ptr;
  typedef map<Symmetry, _ptr, ltsym> _Map;
  typedef typename _Map::value_type _mvtype;
public:
  typedef typename _Map::iterator iterator;
  typedef typename _Map::const_iterator const_iterator;
private:
  int n;
  _Map m;
protected:
  const EquivalenceBundle &ebundle;
public:
  TensorContainer(int nn)
    : n(nn), ebundle(*(tls.ebundle))
  {
  }
  /* This is just a copy constructor. This makes a hard copy of all tensors. */
  TensorContainer(const TensorContainer<_Ttype> &c)
    : n(c.n), m(), ebundle(c.ebundle)
  {
    for (auto it = c.m.begin(); it != c.m.end(); ++it)
      {
        auto *ten = new _Ttype(*((*it).second));
        insert(ten);
      }
  }

  // |TensorContainer| subtensor constructor
  /* This constructor constructs a new tensor container, whose tensors
     are in-place subtensors of the given container. */
  TensorContainer(int first_row, int num, TensorContainer<_Ttype> &c)
    : n(c.n), ebundle(*(tls.ebundle))
  {
    for (auto it = c.m.begin(); it != c.m.end(); ++it)
      {
        auto *t = new _Ttype(first_row, num, *((*it).second));
        insert(t);
      }
  }

  _const_ptr
  get(const Symmetry &s) const
  {
    TL_RAISE_IF(s.num() != num(),
                "Incompatible symmetry lookup in TensorContainer::get");
    auto it = m.find(s);
    if (it == m.end())
      {
        TL_RAISE("Symmetry not found in TensorContainer::get");
        return NULL;
      }
    else
      {
        return (*it).second;
      }
  }

  _ptr
  get(const Symmetry &s)
  {
    TL_RAISE_IF(s.num() != num(),
                "Incompatible symmetry lookup in TensorContainer::get");
    auto it = m.find(s);
    if (it == m.end())
      {
        TL_RAISE("Symmetry not found in TensorContainer::get");
        return NULL;
      }
    else
      {
        return (*it).second;
      }
  }

  bool
  check(const Symmetry &s) const
  {
    TL_RAISE_IF(s.num() != num(),
                "Incompatible symmetry lookup in TensorContainer::check");
    auto it = m.find(s);
    return it != m.end();
  }

  void
  insert(_ptr t)
  {
    TL_RAISE_IF(t->getSym().num() != num(),
                "Incompatible symmetry insertion in TensorContainer::insert");
    TL_RAISE_IF(check(t->getSym()),
                "Tensor already in container in TensorContainer::insert");
    m.insert(_mvtype(t->getSym(), t));
    if (!t->isFinite())
      {
        throw TLException(__FILE__, __LINE__,  "NaN or Inf asserted in TensorContainer::insert");
      }
  }

  void
  remove(const Symmetry &s)
  {
    auto it = m.find(s);
    if (it != m.end())
      {
        _ptr t = (*it).second;
        m.erase(it);
        delete t;
      }
  }

  void
  clear()
  {
    while (!m.empty())
      {
        delete (*(m.begin())).second;
        m.erase(m.begin());
      }
  }

  int
  getMaxDim() const
  {
    int res = -1;
    for (auto run = m.begin(); run != m.end(); ++run)
      {
        int dim = (*run).first.dimen();
        if (dim > res)
          res = dim;
      }
    return res;
  }

  /* Debug print. */
  void
  print() const
  {
    printf("Tensor container: nvars=%d, tensors=%D\n", n, m.size());
    for (const_iterator it = m.begin(); it != m.end(); ++it)
      {
        printf("Symmetry: ");
        (*it).first.print();
        ((*it).second)->print();
      }
  }

  /* Output to the MAT file. */
  void
  writeMat(mat_t *fd, const char *prefix) const
  {
    for (auto it = begin(); it != end(); ++it)
      {
        char lname[100];
        sprintf(lname, "%s_g", prefix);
        const Symmetry &sym = (*it).first;
        for (int i = 0; i < sym.num(); i++)
          {
            char tmp[10];
            sprintf(tmp, "_%d", sym[i]);
            strcat(lname, tmp);
          }
        ConstTwoDMatrix m(*((*it).second));
        m.writeMat(fd, lname);
      }
  }

  /* Output to the Memory Map. */
  void
  writeMMap(map<string, ConstTwoDMatrix> &mm, const string &prefix) const
  {
    ostringstream lname;
    for (const_iterator it = begin(); it != end(); ++it)
      {
        lname.str(prefix);
        lname << "_g";
        const Symmetry &sym = (*it).first;
        for (int i = 0; i < sym.num(); i++)
          lname << "_" << sym[i];
        mm.insert(make_pair(lname.str(), ConstTwoDMatrix(*((*it).second))));
      }
  }

  /* Here we fetch all tensors given by symmetry and equivalence. We go
     through all equivalence classes, calculate implied symmetry, and
     fetch its tensor storing it in the same order to the vector. */

  vector<_const_ptr>
  fetchTensors(const Symmetry &rsym, const Equivalence &e) const
  {
    vector<_const_ptr> res(e.numClasses());
    int i = 0;
    for (auto it = e.begin();
         it != e.end(); ++it, i++)
      {
        Symmetry s(rsym, *it);
        res[i] = get(s);
      }
    return res;
  }

  virtual ~TensorContainer()
  {
    clear();
  }

  int
  num() const
  {
    return n;
  }
  const EquivalenceBundle &
  getEqBundle() const
  {
    return ebundle;
  }

  const_iterator
  begin() const
  {
    return m.begin();
  }
  const_iterator
  end() const
  {
    return m.end();
  }
  iterator
  begin()
  {
    return m.begin();
  }
  iterator
  end()
  {
    return m.end();
  }
};

/* Here is a container storing |UGSTensor|s. We declare |multAndAdd| method. */

class FGSContainer;
class UGSContainer : public TensorContainer<UGSTensor>
{
public:
  UGSContainer(int nn)
    : TensorContainer<UGSTensor>(nn)
  {
  }
  UGSContainer(const UGSContainer &uc)
     
  = default;
  UGSContainer(const FGSContainer &c);
  void multAndAdd(const UGSTensor &t, UGSTensor &out) const;
};

/* Here is a container storing |FGSTensor|s. We declare two versions of
   |multAndAdd| method. The first works for folded $B$ and folded $h$
   tensors, the second works for folded $B$ and unfolded $h$. There is no
   point to do it for unfolded $B$ since the algorithm go through all the
   indices of $B$ and calculates corresponding columns. So, if $B$ is
   needed unfolded, it is more effective to calculate its folded version
   and then unfold by conversion.

   The static member |num_one_time| is a number of columns formed from
   product of $g$ tensors at one time. This is subject to change, probably
   we will have to do some tuning and decide about this number based on
   symmetries, and dimensions in the runtime. */

class FGSContainer : public TensorContainer<FGSTensor>
{
  static const int num_one_time;
public:
  FGSContainer(int nn)
    : TensorContainer<FGSTensor>(nn)
  {
  }
  FGSContainer(const FGSContainer &fc)
     
  = default;
  FGSContainer(const UGSContainer &c);
  void multAndAdd(const FGSTensor &t, FGSTensor &out) const;
  void multAndAdd(const UGSTensor &t, FGSTensor &out) const;
private:
  static Tensor::index getIndices(int num, vector<IntSequence> &out,
                                  const Tensor::index &start,
                                  const Tensor::index &end);
};

#endif
