/*
 * Copyright © 2004 Ondra Kamenik
 * Copyright © 2019-2023 Dynare Team
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

// Tensor containers.

/* One of primary purposes of the tensor library is to perform one step
   of the Faà Di Bruno formula:

                                         ₗ
    [B_sᵏ]_α₁…αₗ = [h_zˡ]_γ₁…γₗ    ∑     ∏  [g_{s^|cₘ|}]_cₘ(α)^γₘ
                                c∈ℳₗ,ₖ ᵐ⁼¹

   where [h_zˡ] and [g_sⁱ] are tensors, ℳₗ,ₖ is the set of all equivalences
   with l classes of k element set, cₘ is an m-class of equivalence c, and |cₘ|
   is its cardinality. Further, cₘ(α) is the sequence of α’s picked by
   equivalence class cₘ.

   In order to accomplish this operation, we basically need some storage of all
   tensors of the form [g_sⁱ]. Note that s can be compound, for instance
   s=(y,u). Then we need storage for [g_y³], [g_y²u], [g_yu⁵], etc.

   We need an object holding all tensors of the same type. Here type means an
   information, that coordinates of the tensors can be of type y or u. We will
   group only tensors, whose symmetry is described by the Symmetry class. These
   are only y²u³, not yuyu². So, we are going to define a class which will hold
   tensors whose symmetries are of type Symmetry and have the same symmetry
   length (number of different coordinate types). Also, for each symmetry there
   will be at most one tensor.

   The class has two purposes. The first is to provide storage (insert and
   retrieve). The second is to perform the above step of Faà Di Bruno. This is
   going through all equivalences with $l$ classes, perform the tensor product
   and add to the result.

   We define a template class TensorContainer. From different instantiations of
   the template class we will inherit to create concrete classes, for example
   container of unfolded general symmetric tensors. The one step of the Faà Di
   Bruno (we call it multAndAdd()) is implemented in the concrete subclasses,
   because the implementation depends on storage. Note even, that multAndAdd()
   has not a template common declaration. This is because sparse tensor h is
   multiplied by folded tensors g yielding folded tensor B, but unfolded tensor
   h is multiplied by unfolded tensors g yielding unfolded tensor B. */

#ifndef T_CONTAINER_H
#define T_CONTAINER_H

#include "Vector.hh"
#include "equivalence.hh"
#include "gs_tensor.hh"
#include "rfs_tensor.hh"
#include "sparse_tensor.hh"
#include "symmetry.hh"
#include "tl_exception.hh"
#include "tl_static.hh"

#include <map>
#include <memory>
#include <string>
#include <utility>

// ltsym predicate
/* We need a predicate on strict weak ordering of
   symmetries. */
struct ltsym
{
  bool
  operator()(const Symmetry& s1, const Symmetry& s2) const
  {
    return s1 < s2;
  }
};

/* Here we define the template class for tensor container. We implement it as
   an stl::map. It is a unique container, no two tensors with same symmetries
   can coexist. Keys of the map are symmetries, values are pointers to tensor.
   The class is responsible for deallocating all tensors. Creation of the
   tensors is done outside.

   The class has integer ‘n’ as its member. It is a number of different
   coordinate types of all contained tensors. Besides intuitive insert and
   retrieve interface, we define a method fetchTensors(), which for a given
   symmetry and given equivalence calculates symmetries implied by the symmetry
   and all equivalence classes, and fetches corresponding tensors in a vector.

   Also, each instance of the container has a reference to EquivalenceBundle
   which allows an access to equivalences. */

template<class _Ttype>
class TensorContainer
{
protected:
  using _Map = std::map<Symmetry, std::unique_ptr<_Ttype>, ltsym>;

private:
  int n;
  _Map m;

public:
  TensorContainer(int nn) : n(nn)
  {
  }
  /* This is just a copy constructor. This makes a hard copy of all tensors. */
  TensorContainer(const TensorContainer<_Ttype>& c) : n(c.n)
  {
    for (const auto& it : c.m)
      insert(std::make_unique<_Ttype>(*(it.second)));
  }
  TensorContainer(TensorContainer<_Ttype>&&) = default;

  // TensorContainer subtensor constructor
  /* This constructor constructs a new tensor container, whose tensors
     are in-place subtensors of the given container. */
  TensorContainer(int first_row, int num, TensorContainer<_Ttype>& c) : n(c.n)
  {
    for (const auto& it : c.m)
      insert(std::make_unique<_Ttype>(first_row, num, *(it.second)));
  }

  TensorContainer<_Ttype>&
  operator=(const TensorContainer<_Ttype>& c)
  {
    n = c.n;
    m.clear();
    for (const auto& it : c.m)
      insert(std::make_unique<_Ttype>(*(it.second)));
  }
  TensorContainer<_Ttype>& operator=(TensorContainer<_Ttype>&&) = default;

  [[nodiscard]] const _Ttype&
  get(const Symmetry& s) const
  {
    TL_RAISE_IF(s.num() != num(), "Incompatible symmetry lookup in TensorContainer::get");
    auto it = m.find(s);
    TL_RAISE_IF(it == m.end(), "Symmetry not found in TensorContainer::get");
    return *(it->second);
  }

  _Ttype&
  get(const Symmetry& s)
  {
    TL_RAISE_IF(s.num() != num(), "Incompatible symmetry lookup in TensorContainer::get");
    auto it = m.find(s);
    TL_RAISE_IF(it == m.end(), "Symmetry not found in TensorContainer::get");
    return *(it->second);
  }

  [[nodiscard]] bool
  check(const Symmetry& s) const
  {
    TL_RAISE_IF(s.num() != num(), "Incompatible symmetry lookup in TensorContainer::check");
    auto it = m.find(s);
    return it != m.end();
  }

  virtual void
  insert(std::unique_ptr<_Ttype> t)
  {
    TL_RAISE_IF(t->getSym().num() != num(),
                "Incompatible symmetry insertion in TensorContainer::insert");
    TL_RAISE_IF(check(t->getSym()), "Tensor already in container in TensorContainer::insert");
    if (!t->isFinite())
      throw TLException(__FILE__, __LINE__, "NaN or Inf asserted in TensorContainer::insert");
    m.emplace(t->getSym(), std::move(t));
  }

  void
  remove(const Symmetry& s)
  {
    m.erase(s);
  }

  void
  clear()
  {
    m.clear();
  }

  [[nodiscard]] int
  getMaxDim() const
  {
    int res = -1;
    for (const auto& run : m)
      if (int dim {run.first.dimen()}; dim > res)
        res = dim;
    return res;
  }

  /* Debug print. */
  void
  print() const
  {
    std::cout << "Tensor container: nvars=" << n << ", tensors=" << m.size() << '\n';
    for (auto& it : *this)
      {
        std::cout << "Symmetry: ";
        it.first.print();
        it.second->print();
      }
  }

  /* Output to the Memory Map. */
  void
  writeMMap(std::map<std::string, ConstTwoDMatrix>& mm, const std::string& prefix) const
  {
    for (auto& it : *this)
      {
        std::string lname = prefix + "_g";
        const Symmetry& sym = it.first;
        for (int i = 0; i < sym.num(); i++)
          lname += '_' + std::to_string(sym[i]);
        mm.emplace(lname, ConstTwoDMatrix(*(it.second)));
      }
  }

  /* Here we fetch all tensors given by symmetry and equivalence. We go
     through all equivalence classes, calculate implied symmetry, and
     fetch its tensor storing it in the same order to the vector. */

  [[nodiscard]] std::vector<const _Ttype*>
  fetchTensors(const Symmetry& rsym, const Equivalence& e) const
  {
    std::vector<const _Ttype*> res(e.numClasses());
    int i = 0;
    for (auto it = e.begin(); it != e.end(); ++it, i++)
      {
        Symmetry s(rsym, *it);
        res[i] = &get(s);
      }
    return res;
  }

  virtual ~TensorContainer() = default;

  [[nodiscard]] int
  num() const
  {
    return n;
  }

  [[nodiscard]] auto
  begin() const
  {
    return m.begin();
  }
  [[nodiscard]] auto
  end() const
  {
    return m.end();
  }
  auto
  begin()
  {
    return m.begin();
  }
  auto
  end()
  {
    return m.end();
  }
};

/* Here is a container storing UGSTensor’s. We declare multAndAdd() method. */

class FGSContainer;
class UGSContainer : public TensorContainer<UGSTensor>
{
public:
  UGSContainer(int nn) : TensorContainer<UGSTensor>(nn)
  {
  }
  UGSContainer(const FGSContainer& c);
  void multAndAdd(const UGSTensor& t, UGSTensor& out) const;
};

/* Here is a container storing FGSTensor’s. We declare two versions of
   multAndAdd() method. The first works for folded B and folded h tensors, the
   second works for folded B and unfolded h. There is no point to do it for
   unfolded B since the algorithm goes through all the indices of B and
   calculates corresponding columns. So, if B is needed unfolded, it is more
   effective to calculate its folded version and then unfold by conversion.

   The static member ‘num_one_time’ is a number of columns formed from product
   of g tensors at one time. This is subject to change, probably we will have
   to do some tuning and decide about this number based on symmetries, and
   dimensions in the runtime. */

class FGSContainer : public TensorContainer<FGSTensor>
{
  static constexpr int num_one_time = 10;

public:
  FGSContainer(int nn) : TensorContainer<FGSTensor>(nn)
  {
  }
  FGSContainer(const UGSContainer& c);
  void multAndAdd(const FGSTensor& t, FGSTensor& out) const;
  void multAndAdd(const UGSTensor& t, FGSTensor& out) const;

private:
  static Tensor::index getIndices(int num, std::vector<IntSequence>& out,
                                  const Tensor::index& start, const Tensor::index& end);
};

#endif
