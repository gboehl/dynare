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
 * along with Dynare.  If not, see <http://www.gnu.org/licenses/>.
 */

#include "kord_exception.hh"
#include "first_order.hh"

#include <dynlapack.h>

double FirstOrder::qz_criterium_global;
std::mutex FirstOrder::mut;

/* This is a function which selects the eigenvalues pair used by LAPACK’s
   dgges. Here we want to select the pairs for which α<β (up to the QZ
   criterium). */

lapack_int
FirstOrder::order_eigs(const double *alphar, const double *alphai, const double *beta)
{
  return (*alphar **alphar + *alphai **alphai < *beta **beta * qz_criterium_global * qz_criterium_global);
}

/* Here we solve the linear approximation. The result are the matrices
   g_y* and gᵤ. The method solves the first derivatives of g so
   that the following equation would be true:

    𝔼ₜ[F(y*ₜ₋₁,uₜ,uₜ₊₁,σ)] = 𝔼ₜ[f(g**(g*(y*ₜ₋₁,uₜ,σ), uₜ₊₁, σ), g(y*ₜ₋₁,uₜ,σ), y*ₜ₋₁,uₜ)]=0

   where f is a given system of equations.

   It is known that g_y* is given by F_y*=0, gᵤ is given by Fᵤ=0, and g_σ is
   zero. The only input to the method are the derivatives ‘fd’ of the system f,
   and partitioning of the vector y (from object). */

void
FirstOrder::solve(const TwoDMatrix &fd)
{
  JournalRecordPair pa(journal);
  pa << "Recovering first order derivatives " << endrec;

  // Solve derivatives ‘gy’

  /* The derivatives g_y* are retrieved from the equation F_y*=0. The
     calculation proceeds as follows:

     1. For each variable appearing at both t-1 and t+1 we add a dummy
        variable, so that the predetermined variables and forward looking would
        be disjoint. This is, the matrix of the first derivatives of the
        system written as:

         [ f_y**₊  f_ys  f_yp  f_yb  f_yf  f_y*₋ ]

        where f_ys, f_yp, f_yb, and f_yf are derivatives w.r.t. static,
        predetermined, both, forward looking at time t, is rewritten as:

         ⎡ f_y**₊  f_ys  f_yp  f_yb   0  f_yf  f_y*₋ ⎤
         ⎣    0      0     0     I   −I    0     0   ⎦

        where the second line has number of rows equal to the number of both
        variables.

     2. Next, provided that forward looking and predetermined are disjoint, the
        equation F_y*=0 is written as:

         [f_y**₊][g**_y*][g*_y*] + [f_ys][gˢ_y*] + [f_y*][g*_y*] + [f_y**][g**_y*] + [f_y*₋] = 0

        This is rewritten as

                        ⎡   I  ⎤                            ⎡   I  ⎤
         [f_y* 0 f_y**₊]⎢ gˢ_y*⎥[g*_y*] + [f_y*₋ f_ys f_y**]⎢ gˢ_y*⎥ = 0
                        ⎣g**_y*⎦                            ⎣g**_y*⎦

        Now, in the above equation, there are the auxiliary variables standing
        for copies of both variables at time t+1. This equation is then
        rewritten as:

                              ⎡   I  ⎤                              ⎡   I  ⎤
         ⎡f_yp f_yb  0 f_y**₊⎤⎢ gˢ_y*⎥[g*_y*] + ⎡f_y*₋ f_ys  0 f_yf⎤⎢ gˢ_y*⎥ = 0
         ⎣  0    I   0    0  ⎦⎣g**_y*⎦          ⎣  0     0  −I   0 ⎦⎣g**_y*⎦

        The two matrices are denoted as D and −E, so the equation takes the form:

           ⎡   I  ⎤            ⎡   I  ⎤
         D ⎢ gˢ_y*⎥[g*_y*] = E ⎢ gˢ_y*⎥
           ⎣g**_y*⎦            ⎣g**_y*⎦

     3. Next we solve the equation by Generalized Schur decomposition:

        ⎡T₁₁ T₁₂⎤⎡Z₁₁ᵀ Z₂₁ᵀ⎤⎡I⎤          ⎡S₁₁ S₁₂⎤⎡Z₁₁ᵀ Z₂₁ᵀ⎤⎡I⎤
        ⎣ 0  T₂₂⎦⎣Z₁₂ᵀ Z₂₂ᵀ⎦⎣X⎦[g*_y*] = ⎣ 0  S₂₂⎦⎣Z₁₂ᵀ Z₂₂ᵀ⎦⎣X⎦

        We reorder the eigenvalue pair so that Sᵢᵢ/Tᵢᵢ with modulus less than
        one would be in the left-upper part.

     4. The Blanchard-Kahn stability argument implies that the pairs with
        modulus less that one will be in and only in S₁₁/T₁₁. The exploding
        paths will be then eliminated when

         ⎡Z₁₁ᵀ Z₂₁ᵀ⎤⎡I⎤   ⎡Y⎤
         ⎣Z₁₂ᵀ Z₂₂ᵀ⎦⎣X⎦ = ⎣0⎦

        From this we have, Y=Z₁₁⁻¹, and X=Z₂₁Y, or equivalently X=−Z₂₂⁻ᵀZ₁₂ᵀ.
        From the equation, we get [g*_y*]=Y⁻¹T₁₁⁻¹S₁₁Y, which is
        Z₁₁T₁₁⁻¹S₁₁Z₁₁⁻¹.

     5. We then copy the derivatives to storage ‘gy’. Note that the derivatives
        of both variables are in X and in [g*_y*], so we check whether the two
        submatrices are the same. The difference is only numerical error.

  */

  // Setup submatrices of ‘f’
  /* Here we setup submatrices of the derivatives ‘fd’. */
  int off = 0;
  ConstTwoDMatrix fyplus(fd, off, ypart.nyss());
  off += ypart.nyss();
  ConstTwoDMatrix fyszero(fd, off, ypart.nstat);
  off += ypart.nstat;
  ConstTwoDMatrix fypzero(fd, off, ypart.npred);
  off += ypart.npred;
  ConstTwoDMatrix fybzero(fd, off, ypart.nboth);
  off += ypart.nboth;
  ConstTwoDMatrix fyfzero(fd, off, ypart.nforw);
  off += ypart.nforw;
  ConstTwoDMatrix fymins(fd, off, ypart.nys());
  off += ypart.nys();
  ConstTwoDMatrix fuzero(fd, off, nu);
  off += nu;

  // Form matrix D
  lapack_int n = ypart.ny()+ypart.nboth;
  TwoDMatrix matD(n, n);
  matD.zeros();
  matD.place(fypzero, 0, 0);
  matD.place(fybzero, 0, ypart.npred);
  matD.place(fyplus, 0, ypart.nys()+ypart.nstat);
  for (int i = 0; i < ypart.nboth; i++)
    matD.get(ypart.ny()+i, ypart.npred+i) = 1.0;
  lapack_int ldb = matD.getLD();

  // Form matrix E
  TwoDMatrix matE(n, n);
  matE.zeros();
  matE.place(fymins, 0, 0);
  matE.place(fyszero, 0, ypart.nys());
  matE.place(fyfzero, 0, ypart.nys()+ypart.nstat+ypart.nboth);
  for (int i = 0; i < ypart.nboth; i++)
    matE.get(ypart.ny()+i, ypart.nys()+ypart.nstat+i) = -1.0;
  matE.mult(-1.0);
  lapack_int lda = matE.getLD();

  // Solve generalized Schur decomposition
  TwoDMatrix vsl(n, n);
  TwoDMatrix vsr(n, n);
  lapack_int ldvsl = vsl.getLD(), ldvsr = vsr.getLD();
  lapack_int lwork = 100*n+16;
  Vector work(lwork);
  auto bwork = std::make_unique<lapack_int[]>(n);
  lapack_int info;
  lapack_int sdim2 = sdim;
  {
    std::lock_guard<std::mutex> lk{mut};
    qz_criterium_global = qz_criterium;
    dgges("N", "V", "S", order_eigs, &n, matE.getData().base(), &lda,
          matD.getData().base(), &ldb, &sdim2, alphar.base(), alphai.base(),
          beta.base(), vsl.getData().base(), &ldvsl, vsr.getData().base(), &ldvsr,
          work.base(), &lwork, bwork.get(), &info);
  }
  if (info)
    throw KordException(__FILE__, __LINE__,
                        "DGGES returns an error in FirstOrder::solve");
  sdim = sdim2;
  bk_cond = (sdim == ypart.nys());

  // Setup submatrices of Z
  ConstGeneralMatrix z11(vsr, 0, 0, ypart.nys(), ypart.nys());
  ConstGeneralMatrix z12(vsr, 0, ypart.nys(), ypart.nys(), n-ypart.nys());
  ConstGeneralMatrix z21(vsr, ypart.nys(), 0, n-ypart.nys(), ypart.nys());
  ConstGeneralMatrix z22(vsr, ypart.nys(), ypart.nys(), n-ypart.nys(), n-ypart.nys());

  // Calculate derivatives of static and forward
  /* Here we calculate X=−Z₂₂⁻ᵀZ₁₂ᵀ, where X is ‘sfder’ in the code. */
  GeneralMatrix sfder(transpose(z12));
  z22.multInvLeftTrans(sfder);
  sfder.mult(-1);

  // Calculate derivatives of predetermined
  /* Here we calculate g*_y* = Z₁₁T₁₁⁻¹S₁₁Z₁₁⁻¹ = Z₁₁T₁₁⁻¹(Z₁₁⁻ᵀS₁₁ᵀ)ᵀ. */
  ConstGeneralMatrix s11(matE, 0, 0, ypart.nys(), ypart.nys());
  ConstGeneralMatrix t11(matD, 0, 0, ypart.nys(), ypart.nys());
  GeneralMatrix dumm(transpose(s11));
  z11.multInvLeftTrans(dumm);
  GeneralMatrix preder(transpose(dumm));
  t11.multInvLeft(preder);
  preder.multLeft(z11);

  // Copy derivatives to ‘gy’
  gy.place(preder, ypart.nstat, 0);
  GeneralMatrix sder(sfder, 0, 0, ypart.nstat, ypart.nys());
  gy.place(sder, 0, 0);
  GeneralMatrix fder(sfder, ypart.nstat+ypart.nboth, 0, ypart.nforw, ypart.nys());
  gy.place(fder, ypart.nstat+ypart.nys(), 0);

  // Check difference for derivatives of both
  GeneralMatrix bder(const_cast<const GeneralMatrix &>(sfder), ypart.nstat, 0, ypart.nboth, ypart.nys());
  GeneralMatrix bder2(preder, ypart.npred, 0, ypart.nboth, ypart.nys());
  bder.add(-1, bder2);
  b_error = bder.getData().getMax();

  // Solve derivatives ‘gu’

  /* The equation Fᵤ=0 can be written as

      [f_y**₊][g**_y*][gᵤ*] + [f_y][gᵤ] + [fᵤ] = 0

     and rewritten as

      [f_y + [0 f_y**₊·g**_y* 0] ] gᵤ = fᵤ

     This is exactly what is done here. The matrix [f_y + [0 f_y**₊·g**_y* 0] ]
     is ‘matA’ in the code. */
  GeneralMatrix matA(ypart.ny(), ypart.ny());
  matA.zeros();
  ConstGeneralMatrix gss(gy, ypart.nstat+ypart.npred, 0, ypart.nyss(), ypart.nys());
  matA.place(fyplus * gss, 0, ypart.nstat);
  ConstGeneralMatrix fyzero(fd, 0, ypart.nyss(), ypart.ny(), ypart.ny());
  matA.add(1.0, fyzero);
  gu.zeros();
  gu.add(-1.0, fuzero);
  ConstGeneralMatrix(matA).multInvLeft(gu);

  journalEigs();

  KORD_RAISE_IF_X(!bk_cond,
                  "The model is not Blanchard-Kahn stable",
                  KORD_MD_NOT_STABLE);
  if (!gy.isFinite() || !gu.isFinite())
    throw KordException(__FILE__, __LINE__,
                        "NaN or Inf asserted in first order derivatives in FirstOrder::solve()");
}

void
FirstOrder::journalEigs()
{
  if (bk_cond)
    {
      JournalRecord jr(journal);
      jr << "Blanchard-Kahn conditition satisfied, model stable" << endrec;
    }
  else
    {
      JournalRecord jr(journal);
      jr << "Blanchard-Kahn condition not satisfied, model not stable: sdim=" << sdim
         << " npred=" << ypart.nys() << endrec;
    }
  if (!bk_cond)
    for (int i = 0; i < alphar.length(); i++)
      {
        if (i == sdim || i == ypart.nys())
          {
            JournalRecord jr(journal);
            jr << u8"──────────────────────────────────────────────────── ";
            if (i == sdim)
              jr << "sdim";
            else
              jr << "npred";
            jr << endrec;
          }
        JournalRecord jr(journal);
        double mod = std::sqrt(alphar[i]*alphar[i]+alphai[i]*alphai[i]);
        mod = mod/std::round(100000*std::abs(beta[i]))*100000;
        jr << i << "\t(" << alphar[i] << "," << alphai[i] << ") / " << beta[i]
           << "  \t" << mod << endrec;
      }
}
