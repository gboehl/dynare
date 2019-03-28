/* $Header: /var/lib/cvs/dynare_cpp/sylv/cc/TriangularSylvester.h,v 1.1.1.1 2004/06/04 13:01:03 kamenik Exp $ */

/* Tag $Name:  $ */

#ifndef TRIANGULAR_SYLVESTER_H
#define TRIANGULAR_SYLVESTER_H

#include "SylvesterSolver.hh"
#include "KronVector.hh"
#include "QuasiTriangular.hh"
#include "QuasiTriangularZero.hh"
#include "SimilarityDecomp.hh"

#include <memory>

class TriangularSylvester : public SylvesterSolver
{
  const std::unique_ptr<const QuasiTriangular> matrixKK;
  const std::unique_ptr<const QuasiTriangular> matrixFF;
public:
  TriangularSylvester(const QuasiTriangular &k, const QuasiTriangular &f);
  TriangularSylvester(const SchurDecompZero &kdecomp, const SchurDecomp &fdecomp);
  TriangularSylvester(const SchurDecompZero &kdecomp, const SimilarityDecomp &fdecomp);
  
  ~TriangularSylvester() override = default;
  void print() const;
  void solve(SylvParams &pars, KronVector &d) const override;

  void solvi(double r, KronVector &d, double &eig_min) const;
  void solvii(double alpha, double beta1, double beta2,
              KronVector &d1, KronVector &d2,
              double &eig_min) const;
  void solviip(double alpha, double betas,
               KronVector &d, double &eig_min) const;
  /* Computes:
     ⎛x₁⎞ ⎛d₁⎞ ⎛ α −β₁⎞           ⎛d₁⎞
     ⎢  ⎥=⎢  ⎥+⎢      ⎥⊗Fᵀ⊗Fᵀ⊗…⊗K·⎢  ⎥
     ⎝x₂⎠ ⎝d₂⎠ ⎝−β₂ α ⎠           ⎝d₂⎠
  */
  void linEval(double alpha, double beta1, double beta2,
               KronVector &x1, KronVector &x2,
               const ConstKronVector &d1, const ConstKronVector &d2) const;
  void
  linEval(double alpha, double beta1, double beta2,
          KronVector &x1, KronVector &x2,
          const KronVector &d1, const KronVector &d2) const
  {
    linEval(alpha, beta1, beta2, x1, x2,
            ConstKronVector(d1), ConstKronVector(d2));
  }

  /* Computes:
     ⎛x₁⎞ ⎛d₁⎞   ⎛γ −δ₁⎞           ⎛d₁⎞       ⎛γ −δ₁⎞²              ⎛d₁⎞
     ⎢  ⎥=⎢  ⎥+2α⎢     ⎥⊗Fᵀ⊗Fᵀ⊗…⊗K·⎢  ⎥+(α²+β)⎢     ⎥ ⊗Fᵀ²⊗Fᵀ²⊗…⊗K²·⎢  ⎥
     ⎝x₂⎠ ⎝d₂⎠   ⎝δ₂ γ ⎠           ⎝d₂⎠       ⎝δ₂ γ ⎠               ⎝d₂⎠
  */
  void quaEval(double alpha, double betas,
               double gamma, double delta1, double delta2,
               KronVector &x1, KronVector &x2,
               const ConstKronVector &d1, const ConstKronVector &d2) const;
  void
  quaEval(double alpha, double betas,
          double gamma, double delta1, double delta2,
          KronVector &x1, KronVector &x2,
          const KronVector &d1, const KronVector &d2) const
  {
    quaEval(alpha, betas, gamma, delta1, delta2, x1, x2,
            ConstKronVector(d1), ConstKronVector(d2));
  }
private:
  /* Returns square of size of minimal eigenvalue of the system solved,
     now obsolete */
  double getEigSep(int depth) const;
  // Recursively calculates kronecker product of complex vectors (used in getEigSep)
  static void multEigVector(KronVector &eig, const Vector &feig, const Vector &keig);

  using const_diag_iter = QuasiTriangular::const_diag_iter;
  using const_row_iter = QuasiTriangular::const_row_iter;

  // Called from solvi
  void solviRealAndEliminate(double r, const_diag_iter di,
                             KronVector &d, double &eig_min) const;
  void solviComplexAndEliminate(double r, const_diag_iter di,
                                KronVector &d, double &eig_min) const;
  // Called from solviip
  void solviipRealAndEliminate(double alpha, double betas,
                               const_diag_iter di, const_diag_iter dsi,
                               KronVector &d, double &eig_min) const;
  void solviipComplexAndEliminate(double alpha, double betas,
                                  const_diag_iter di, const_diag_iter dsi,
                                  KronVector &d, double &eig_min) const;
  // Eliminations
  void solviEliminateReal(const_diag_iter di, KronVector &d,
                          const KronVector &y, double divisor) const;
  void solviEliminateComplex(const_diag_iter di, KronVector &d,
                             const KronVector &y1, const KronVector &y2,
                             double divisor) const;
  void solviipEliminateReal(const_diag_iter di, const_diag_iter dsi,
                            KronVector &d,
                            const KronVector &y1, const KronVector &y2,
                            double divisor, double divisor2) const;
  void solviipEliminateComplex(const_diag_iter di, const_diag_iter dsi,
                               KronVector &d,
                               const KronVector &y1, const KronVector &y11,
                               const KronVector &y2, const KronVector &y22,
                               double divisor) const;
  // Lemma 2
  void solviipComplex(double alpha, double betas, double gamma,
                      double delta1, double delta2,
                      KronVector &d1, KronVector &d2,
                      double &eig_min) const;
  // Norms for what we consider zero on diagonal of F
  static constexpr double diag_zero = 1.e-15;
  static constexpr double diag_zero_sq = diag_zero*diag_zero;
};

#endif /* TRIANGULAR_SYLVESTER_H */
