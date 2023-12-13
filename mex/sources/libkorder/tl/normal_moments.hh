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
 * along with Dynare.  If not, see <https://www.gnu.org/licenses/>.
 */

// Moments of normal distribution.

/* Here we calculate the higher order moments of normally distributed
   random vector u with means equal to zero and given
   variance-covariance matrix V, i.e. u↝𝒩(0,V). The moment
   generating function for such distribution is f(t)=e^½tᵀVt. If
   we derivate it w.r.t. t and unfold the higher dimensional tensors
   row-wise, we obtain terms like:

    ∂f(t)/∂t = f(t)·Vt
    ∂²f(t)/∂t² = f(t)·(Vt⊗v)
    ∂³f(t)/∂t³ = f(t)·(Vt⊗Vt⊗Vt + P_?(v⊗Vt) + P_?(Vt⊗v) + v⊗Vt)
    ∂⁴f(t)/∂t⁴ = f(t)·(Vt⊗Vt⊗Vt⊗Vt + S_?(v⊗Vt⊗Vt) + S_?(Vt⊗v⊗Vt) + S_?(Vt⊗Vt⊗v) + S_?(v⊗v))

   where v is vectorized V (v=vec(V)), and P_? is a suitable row permutation
   (corresponds to permutation of multidimensional indices) which permutes the
   tensor data, so that the index of a variable being derived would be the
   last. This ensures that all (permuted) tensors can be summed yielding a
   tensor whose indices have some order (in here we chose the order that more
   recent derivating variables are to the right). Finally, S_? is a suitable sum
   of various P_?.

   We are interested in S_? multiplying the Kronecker powers ⊗ⁿv. The S_? is a
   (possibly) multi-set of permutations of even order. Note that we know the
   number of permutations in S_?. The above formulas for f(t) derivatives are
   valid also for monomial u, and from literature we know that 2n-th moment is
   (2n!)/(n!2ⁿ)·σ². So there are (2n!)/(n!2ⁿ) permutations in S_?.

   In order to find the S_? permutation we need to define a couple of things.
   First we define a sort of equivalence between the permutations applicable to
   even number of indices. We write P₁ ≡ P₂ whenever P₁⁻¹∘P₂ permutes only whole
   pairs, or items within pairs, but not indices across the pairs. For instance
   the permutations (0,1,2,3) and (3,2,0,1) are equivalent, but (0,2,1,3) is
   not equivalent with the two. Clearly, the ≡ relationship is an equivalence.

   This allows to define a relation ⊑ between the permutation multi-sets S,
   which is basically the subset relation ⊆ but with respect to the equivalence
   relation ≡, more formally:

    S₁ ⊑ S₂  iff  P ∈ S₁ ⇒ ∃Q ∈ S₂ : P ≡ Q

   This induces an equivalence S₁ ≡ S₂.

   Now let Fₙ denote a set of permutations on 2n indices which is maximal with
   respect to ⊑, and minimal with respect to ≡ (in other words, it contains
   everything up to the equivalence ≡). It is straightforward to calculate a
   number of permutations in Fₙ. This is a total number of all permutations
   of 2n divided by permutations of pairs divided by permutations within the
   pairs. This is (2n!)/(n!2ⁿ).

   We now prove that S_? ≡ Fₙ. Clearly S_? ⊑ Fₙ, since Fₙ is maximal. In order
   to prove that Fₙ ⊑ S_?, let us assert that for any permutation P and for
   any (semi-)positive definite matrix V we have
   P·S_?(⊗ⁿv)=S_?(⊗ⁿv). Below we show that there is a positive
   definite matrix V of some dimension that for any two permutation
   multi-sets S₁, S₂, we have

    S₁ ≢ S₂  ⇒  S₁(⊗ⁿv) ≠ S₂(⊗ⁿv)

   So it follows that for any permutation P, we have P·S_? ≡ S_?. For a purpose
   of contradiction let P ∈ Fₙ be a permutation which is not equivalent to any
   permutation from S_?. Since S_? is non-empty, let us pick P₀ ∈ S_?. Now
   assert that P₀⁻¹S_? ≢ P⁻¹S_? since the first contains an identity and the
   second does not contain a permutation equivalent to identity. Thus we have
   (P∘P₀⁻¹)S_? ≢ S_? which gives the contradiction and we have proved that
   Fₙ ⊑ S_?. Thus Fₙ ≡ S_?. Moreover, we know that S_? and Fₙ have the same
   number of permutations, hence the minimality of S_? with respect to ≡.

   Now it suffices to prove that there exists a positive definite V such that
   for any two permutation multi-sets S₁ and S₂ holds S₁ ≢ S₂ ⇒ S₁(⊗ⁿv) ≠ S₂⊗ⁿv.
   If V is a n×n matrix, then S₁ ≢ S₂ implies that there is
   identically nonzero polynomial of elements from V of order n over
   integers. If V=AᵀA then there is identically non-zero polynomial of
   elements from A of order 2n. This means, that we have to find n(n+1)/2
   tuple x of real numbers such that all identically non-zero polynomials p
   of order 2n over integers yield p(x)≠0.

   The x is constructed as follows: x_i = π^log(rᵢ), where rᵢ is i-th prime.
   Let us consider the monom x₁^j₁·…·xₖ^jₖ. When the monom is evaluated, we get

    π^{log(r₁^j₁)+…+log(rₖ^jₖ)}= π^log(r₁^j₁+…+rₖ^jₖ)

   Now it is easy to see that if an integer combination of such terms is zero,
   then the combination must be either trivial or sum to 0 and all monoms
   must be equal. Both cases imply a polynomial identically equal to zero. So,
   any non-trivial integer polynomial evaluated at x must be non-zero.

   So, having this result in hand, now it is straightforward to calculate
   higher moments of normal distribution. Here we define a container, which
   does the job. In its constructor, we simply calculate Kronecker powers of v
   and apply Fₙ to ⊗ⁿv. Fₙ is, in fact, a set of all equivalences in sense of
   class Equivalence over 2n elements, having n classes each of them having
   exactly 2 elements.
*/

#ifndef NORMAL_MOMENTS_HH
#define NORMAL_MOMENTS_HH

#include "t_container.hh"

class UNormalMoments : public TensorContainer<URSingleTensor>
{
public:
  UNormalMoments(int maxdim, const TwoDMatrix& v);

private:
  void generateMoments(int maxdim, const TwoDMatrix& v);
  static bool selectEquiv(const Equivalence& e);
};

class FNormalMoments : public TensorContainer<FRSingleTensor>
{
public:
  FNormalMoments(const UNormalMoments& moms);
};

#endif
