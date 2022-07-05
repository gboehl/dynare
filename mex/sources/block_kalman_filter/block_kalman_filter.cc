/*
 * Copyright © 2007-2022 Dynare Team
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

#include <cmath>
#include <algorithm>
#include <omp.h>

#include "block_kalman_filter.hh"

#define BLAS
//#define CUBLAS

#ifdef CUBLAS
# include <cuda_runtime.h>
# include <cublas_v2.h>
#endif

void
mexDisp(const mxArray *P)
{
  size_t n = mxGetN(P);
  size_t m = mxGetM(P);
  const double *M = mxGetPr(P);
  mexPrintf("%d x %d\n", m, n);
  mexEvalString("drawnow;");
  for (size_t i = 0; i < m; i++)
    {
      for (size_t j = 0; j < n; j++)
        mexPrintf(" %9.4f", M[i+ j * m]);
      mexPrintf("\n");
    }
  mexEvalString("drawnow;");
}

void
mexDisp(const double *M, int m, int n)
{
  mexPrintf("%d x %d\n", m, n);
  mexEvalString("drawnow;");
  for (int i = 0; i < m; i++)
    {
      for (int j = 0; j < n; j++)
        mexPrintf(" %9.4f", M[i+ j * m]);
      mexPrintf("\n");
    }
  mexEvalString("drawnow;");
}

/*if block
  %nz_state_var = M_.nz_state_var;
  while notsteady && t<smpl
  t  = t+1;
  v  = Y(:,t)-a(mf);
  F  = P(mf,mf) + H;
  if rcond(F) < kalman_tol
  if ~all(abs(F(:))<kalman_tol)
  return
  else
  a = T*a;
  P = T*P*transpose(T)+QQ;
  end
  else
  F_singular = 0;
  dF     = det(F);
  iF     = inv(F);
  lik(t) = log(dF)+transpose(v)*iF*v;
  K      = P(:,mf)*iF;
  a      = T*(a+K*v);
  P = block_pred_Vcov_KF(mf, P, K, T, QQ);
  %P      = T*(P-K*P(mf,:))*transpose(T)+QQ;
  notsteady = max(abs(K(:)-oldK)) > riccati_tol;
  oldK = K(:);
  end
  end;
  else
  while notsteady && t<smpl
  t  = t+1;
  v  = Y(:,t)-a(mf);
  F  = P(mf,mf) + H;
  if rcond(F) < kalman_tol
  if ~all(abs(F(:))<kalman_tol)
  return
  else
  a = T*a;
  P = T*P*transpose(T)+QQ;
  end
  else
  F_singular = 0;
  dF     = det(F);
  iF     = inv(F);
  lik(t) = log(dF)+transpose(v)*iF*v;
  K      = P(:,mf)*iF;
  a      = T*(a+K*v);
  P      = T*(P-K*P(mf,:))*transpose(T)+QQ;
  notsteady = max(abs(K(:)-oldK)) > riccati_tol;
  oldK = K(:);

  end
  end
  end
*/

bool
not_all_abs_F_bellow_crit(const double *F, int size, double crit)
{
  int i = 0;
  while (i < size && abs(F[i]) < crit)
    i++;

  if (i < size)
    return false;
  else
    return true;
}

double
det(const double *F, int dim, const lapack_int *ipiv)
{
  double det = 1.0;
  for (int i = 0; i < dim; i++)
    if (ipiv[i] - 1 == i)
      det *= F[i + i * dim];
    else
      det *= -F[i + i * dim];
  return det;
}

BlockKalmanFilter::BlockKalmanFilter(int nrhs, const mxArray *prhs[])
{
  if (nrhs != 13 && nrhs != 16)
    mexErrMsgTxt("block_kalman_filter requires exactly \n  13 input arguments for standard Kalman filter \nor\n  16 input arguments for missing observations Kalman filter.");
  if (nrhs == 16)
    missing_observations = true;
  else
    missing_observations = false;
  if (missing_observations)
    {
      if (!mxIsCell(prhs[0]))
        mexErrMsgTxt("the first input argument of block_missing_observations_kalman_filter must be a Cell Array.");
      pdata_index = prhs[0];
      if (!mxIsDouble(prhs[1]))
        mexErrMsgTxt("the second input argument of block_missing_observations_kalman_filter must be a scalar.");
      number_of_observations = ceil(mxGetScalar(prhs[1]));
      if (!mxIsDouble(prhs[2]))
        mexErrMsgTxt("the third input argument of block_missing_observations_kalman_filter must be a scalar.");
      no_more_missing_observations = ceil(mxGetScalar(prhs[2]));
      pT = mxDuplicateArray(prhs[3]);
      pR = mxDuplicateArray(prhs[4]);
      pQ = mxDuplicateArray(prhs[5]);
      pH = mxDuplicateArray(prhs[6]);
      pP = mxDuplicateArray(prhs[7]);
      pY = mxDuplicateArray(prhs[8]);
      start = mxGetScalar(prhs[9]);
      mfd = mxGetPr(prhs[10]);
      kalman_tol = mxGetScalar(prhs[11]);
      riccati_tol = mxGetScalar(prhs[12]);
      nz_state_var = mxGetPr(prhs[13]);
      n_diag = mxGetScalar(prhs[14]);
      pure_obs = mxGetScalar(prhs[15]);
    }
  else
    {
      no_more_missing_observations = 0;
      pT = mxDuplicateArray(prhs[0]);
      pR = mxDuplicateArray(prhs[1]);
      pQ = mxDuplicateArray(prhs[2]);
      pH = mxDuplicateArray(prhs[3]);
      pP = mxDuplicateArray(prhs[4]);
      pY = mxDuplicateArray(prhs[5]);
      start = mxGetScalar(prhs[6]);
      /*Defining the initials values*/
      n = mxGetN(pT); // Number of state variables.
      pp = mxGetM(pY); // Maximum number of observed variables.
      smpl = mxGetN(pY); // Sample size.          ;
      mfd = mxGetPr(prhs[7]);
      kalman_tol = mxGetScalar(prhs[8]);
      riccati_tol = mxGetScalar(prhs[9]);
      nz_state_var = mxGetPr(prhs[10]);
      n_diag = mxGetScalar(prhs[11]);
      pure_obs = mxGetScalar(prhs[12]);
    }
  T = mxGetPr(pT);
  R = mxGetPr(pR);
  Q = mxGetPr(pQ);
  H = mxGetPr(pH);
  P = mxGetPr(pP);
  Y = mxGetPr(pY);

  n = mxGetN(pT); // Number of state variables.
  pp = mxGetM(pY); // Maximum number of observed variables.
  smpl = mxGetN(pY); // Sample size.          ;
  n_state = n - pure_obs;

  /*mexPrintf("T\n");
    mexDisp(pT);*/

  H_size = mxGetN(pH) * mxGetM(pH);

  n_shocks = mxGetM(pQ);

  if (missing_observations)
    if (mxGetNumberOfElements(pdata_index) != static_cast<unsigned int>(smpl))
      mexErrMsgTxt("the number of element in the cell array passed to block_missing_observation_kalman_filter as first argument has to be equal to the smpl size");

  i_nz_state_var = std::make_unique<int[]>(n);
  for (int i = 0; i < n; i++)
    i_nz_state_var[i] = nz_state_var[i];

  pa = mxCreateDoubleMatrix(n, 1, mxREAL); // State vector.
  a = mxGetPr(pa);
  tmp_a = std::make_unique<double[]>(n);
  dF = 0.0; // det(F).

  p_tmp1 = mxCreateDoubleMatrix(n, n_shocks, mxREAL);
  tmp1 = mxGetPr(p_tmp1);
  t = 0; // Initialization of the time index.
  plik = mxCreateDoubleMatrix(smpl, 1, mxREAL);
  lik = mxGetPr(plik);
  Inf = mxGetInf();
  LIK = 0.0; // Default value of the log likelihood.
  notsteady = true; // Steady state flag.
  F_singular = true;
  v_pp = std::make_unique<double[]>(pp);
  v_n = std::make_unique<double[]>(n);
  mf = std::make_unique<int[]>(pp);
  for (int i = 0; i < pp; i++)
    mf[i] = mfd[i] - 1;

  /*compute QQ = R*Q*transpose(R)*/ // Variance of R times the vector of structural innovations.;
  // tmp = R * Q;
  for (int i = 0; i < n; i++)
    for (int j = 0; j < n_shocks; j++)
      {
        double res = 0.0;
        for (int k = 0; k < n_shocks; k++)
          res += R[i + k * n] * Q[j * n_shocks + k];
        tmp1[i + j * n] = res;
      }

  // QQ = tmp * transpose(R)
  pQQ = mxCreateDoubleMatrix(n, n, mxREAL);
  QQ = mxGetPr(pQQ);
  for (int i = 0; i < n; i++)
    for (int j = i; j < n; j++)
      {
        double res = 0.0;
        for (int k = 0; k < n_shocks; k++)
          res += tmp1[i + k * n] * R[k * n + j];
        QQ[i + j * n] = QQ[j + i * n] = res;
      }
  mxDestroyArray(p_tmp1);

  pv = mxCreateDoubleMatrix(pp, 1, mxREAL);
  v = mxGetPr(pv);
  pF = mxCreateDoubleMatrix(pp, pp, mxREAL);
  F = mxGetPr(pF);
  piF = mxCreateDoubleMatrix(pp, pp, mxREAL);
  iF = mxGetPr(piF);
  lw = pp * 4;
  w = std::make_unique<double[]>(lw);
  iw = std::make_unique<lapack_int[]>(pp);
  ipiv = std::make_unique<lapack_int[]>(pp);
  info = 0;
#if defined(BLAS) || defined(CUBLAS)
  p_tmp = mxCreateDoubleMatrix(n, n, mxREAL);
  tmp = mxGetPr(p_tmp);
  p_P_t_t1 = mxCreateDoubleMatrix(n, n, mxREAL);
  P_t_t1 = mxGetPr(p_P_t_t1);
  pK = mxCreateDoubleMatrix(n, n, mxREAL);
  K = mxGetPr(pK);
  p_K_P = mxCreateDoubleMatrix(n, n, mxREAL);
  K_P = mxGetPr(p_K_P);
  oldK = std::make_unique<double[]>(n * n);
  P_mf = std::make_unique<double[]>(n * n);
  for (int i = 0; i < n  * n; i++)
    oldK[i] = Inf;
#else
  p_tmp = mxCreateDoubleMatrix(n, n_state, mxREAL);
  tmp = mxGetPr(p_tmp);
  p_P_t_t1 = mxCreateDoubleMatrix(n_state, n_state, mxREAL);
  P_t_t1 = mxGetPr(p_P_t_t1);
  pK = mxCreateDoubleMatrix(n, pp, mxREAL);
  K = mxGetPr(pK);
  p_K_P = mxCreateDoubleMatrix(n_state, n_state, mxREAL);
  K_P = mxGetPr(p_K_P);
  oldK = std::make_unique<double[]>(n * pp);
  P_mf = std::make_unique<double[]>(n * pp);
  for (int i = 0; i < n  * pp; i++)
    oldK[i] = Inf;
#endif
}

void
BlockKalmanFilter::block_kalman_filter_ss()
{
  if (t+1 < smpl)
    while (t < smpl)
      {
        //v = Y(:,t)-a(mf);
        for (int i = 0; i < pp; i++)
          v[i] = Y[i + t * pp] - a[mf[i]];

        //a = T*(a+K*v);
        for (int i = pure_obs; i < n; i++)
          {
            double res = 0.0;
            for (int j = 0; j < pp; j++)
              res += K[j  * n + i] * v[j];
            v_n[i] = res + a[i];
          }
        for (int i = 0; i < n; i++)
          {
            double res = 0.0;
            for (int j = pure_obs; j < n; j++)
              res += T[j  * n + i] * v_n[j];
            a[i] = res;
          }

        //lik(t) = transpose(v)*iF*v;
        for (int i = 0; i < pp; i++)
          {
            double res = 0.0;
            for (int j = 0; j < pp; j++)
              res += v[j] * iF[j * pp + i];
            v_pp[i] = res;
          }
        double res = 0.0;
        for (int i = 0; i < pp; i++)
          res += v_pp[i] * v[i];

        lik[t] = (log(dF) + res + pp * log(2.0*M_PI))/2;
        if (t + 1 > start)
          LIK += lik[t];

        t++;
      }
}

bool
BlockKalmanFilter::block_kalman_filter(int nlhs, mxArray *plhs[])
{
  while (notsteady && t < smpl)
    {
      if (missing_observations)
        {
          // retrieve the d_index
          pd_index = mxGetCell(pdata_index, t);
          dd_index = mxGetPr(pd_index);
          size_d_index = mxGetM(pd_index);
          d_index.resize(size_d_index);
          for (int i = 0; i < size_d_index; i++)
            d_index[i] = ceil(dd_index[i]) - 1;

          //v = Y(:,t) - a(mf)
          int i_i = 0;
          //#pragma omp parallel for shared(v, i_i, d_index)
          for (auto i = d_index.begin(); i != d_index.end(); i++)
            {
              //mexPrintf("i_i=%d, omp_get_max_threads()=%d\n",i_i,omp_get_max_threads());
              v[i_i] = Y[*i + t * pp] - a[mf[*i]];
              i_i++;
            }

          //F  = P(mf,mf) + H;
          i_i = 0;
          if (H_size == 1)
            //#pragma omp parallel for shared(iF, F, i_i)
            for (auto i = d_index.begin(); i != d_index.end(); i++, i_i++)
              {
                int j_j = 0;
                for (auto j = d_index.begin(); j != d_index.end(); j++, j_j++)
                  iF[i_i + j_j * size_d_index] = F[i_i + j_j * size_d_index] = P[mf[*i] + mf[*j] * n] + H[0];
              }
          else
            //#pragma omp parallel for shared(iF, F, P, H, mf, i_i)
            for (auto i = d_index.begin(); i != d_index.end(); i++, i_i++)
              {
                int j_j = 0;
                for (auto j = d_index.begin(); j != d_index.end(); j++, j_j++)
                  iF[i_i + j_j * size_d_index] = F[i_i + j_j * size_d_index] = P[mf[*i] + mf[*j] * n] + H[*i + *j * pp];
              }
        }
      else
        {
          size_d_index = pp;

          //v = Y(:,t) - a(mf)
          for (int i = 0; i < pp; i++)
            v[i] = Y[i + t * pp] - a[mf[i]];

          //F  = P(mf,mf) + H;
          if (H_size == 1)
            for (int i = 0; i < pp; i++)
              for (int j = 0; j < pp; j++)
                iF[i + j * pp] = F[i + j * pp] = P[mf[i] + mf[j] * n] + H[0];
          else
            for (int i = 0; i < pp; i++)
              for (int j = 0; j < pp; j++)
                iF[i + j * pp] = F[i + j * pp] = P[mf[i] + mf[j] * n] + H[i + j * pp];
        }

      /* Computes the norm of iF */
      double anorm = dlange("1", &size_d_index, &size_d_index, iF, &size_d_index, w.get());
      //mexPrintf("anorm = %f\n",anorm);

      /* Modifies F in place with a LU decomposition */
      dgetrf(&size_d_index, &size_d_index, iF, &size_d_index, ipiv.get(), &info);
      if (info != 0)
        mexPrintf("dgetrf failure with error %d\n", static_cast<int>(info));

      /* Computes the reciprocal norm */
      dgecon("1", &size_d_index, iF, &size_d_index, &anorm, &rcond, w.get(), iw.get(), &info);
      if (info != 0)
        mexPrintf("dgecon failure with error %d\n", static_cast<int>(info));

      if (rcond < kalman_tol)
        if (not_all_abs_F_bellow_crit(F, size_d_index * size_d_index, kalman_tol)) //~all(abs(F(:))<kalman_tol)
          {
            mexPrintf("error: F singular\n");
            LIK = Inf;
            if (nlhs == 2)
              for (int i = t; i < smpl; i++)
                lik[i] = Inf;
            // info = 0
            return_results_and_clean(nlhs, plhs);
            return false;
          }
        else
          {
            mexPrintf("F singular\n");

            //a = T*a;
            for (int i = 0; i < n; i++)
              {
                double res = 0.0;
                for (int j = pure_obs; j < n; j++)
                  res += T[i + j *n] * a[j];
                tmp_a[i] = res;
              }
            std::copy_n(tmp_a.get(), n, a);

            //P = T*P*transpose(T)+QQ;
            std::fill_n(tmp, 0, n * n_state);

            for (int i = 0; i < n; i++)
              for (int j = pure_obs; j < n; j++)
                {
                  int j1 = j - pure_obs;
                  int j1_n_state = j1 * n_state - pure_obs;
                  for (int k = pure_obs; k < i_nz_state_var[i]; k++)
                    tmp[i + j1 * n] += T[i + k * n] * P[k + j1_n_state];
                }

            std::fill_n(P, 0, n * n);
            int n_n_obs = n * pure_obs;
            for (int i = 0; i < n; i++)
              for (int j = i; j < n; j++)
                for (int k = pure_obs; k < i_nz_state_var[j]; k++)
                  {
                    int k_n = k * n;
                    P[i * n + j] += tmp[i + k_n - n_n_obs] * T[j + k_n];
                  }

            for (int i = 0; i < n; i++)
              {
                for (int j = i; j < n; j++)
                  P[j + i * n] += QQ[j + i * n];
                for (int j = i + 1; j < n; j++)
                  P[i + j * n] = P[j + i * n];
              }
          }
      else
        {
          F_singular = false;

          //dF     = det(F);
          dF = det(iF, size_d_index, ipiv.get());

          //iF     = inv(F);
          //int lwork = 4/*2*/* pp;
          dgetri(&size_d_index, iF, &size_d_index, ipiv.get(), w.get(), &lw, &info);
          if (info != 0)
            mexPrintf("dgetri failure with error %d\n", static_cast<int>(info));

          //lik(t) = log(dF)+transpose(v)*iF*v;
#pragma omp parallel for shared(v_pp)
          for (int i = 0; i < size_d_index; i++)
            {
              double res = 0.0;
              for (int j = 0; j < size_d_index; j++)
                res += v[j] * iF[j  * size_d_index + i];
              v_pp[i] = res;
            }
          double res = 0.0;
          for (int i = 0; i < size_d_index; i++)
            res += v_pp[i] * v[i];

          lik[t] = (log(dF) + res + size_d_index * log(2.0*M_PI))/2;
          if (t + 1 >= start)
            LIK += lik[t];

          if (missing_observations)
            //K      = P(:,mf)*iF;
#pragma omp parallel for shared(P_mf)
            for (int i = 0; i < n; i++)
              {
                int j_j = 0;
                //for (int j = 0; j < pp; j++)
                for (auto j = d_index.begin(); j != d_index.end(); j++, j_j++)
                  P_mf[i + j_j * n] = P[i + mf[*j] * n];
              }
          else
            //K      = P(:,mf)*iF;
            for (int i = 0; i < n; i++)
              for (int j = 0; j < pp; j++)
                P_mf[i + j * n] = P[i + mf[j] * n];

#pragma omp parallel for shared(K)
          for (int i = 0; i < n; i++)
            for (int j = 0; j < size_d_index; j++)
              {
                double res = 0.0;
                int j_pp = j * size_d_index;
                for (int k = 0; k < size_d_index; k++)
                  res += P_mf[i + k * n] * iF[j_pp + k];
                K[i + j * n] = res;
              }

          //a      = T*(a+K*v);
#pragma omp parallel for shared(v_n)
          for (int i = pure_obs; i < n; i++)
            {
              double res = 0.0;
              for (int j = 0; j < size_d_index; j++)
                res += K[j  * n + i] * v[j];
              v_n[i] = res + a[i];
            }

#pragma omp parallel for shared(a)
          for (int i = 0; i < n; i++)
            {
              double res = 0.0;
              for (int j = pure_obs; j < n; j++)
                res += T[j  * n + i] * v_n[j];
              a[i] = res;
            }

          if (missing_observations)
            {
              //P      = T*(P-K*P(mf,:))*transpose(T)+QQ;
              int i_i = 0;
              //#pragma omp parallel for shared(P_mf)
              for (auto i = d_index.begin(); i != d_index.end(); i++, i_i++)
                for (int j = pure_obs; j < n; j++)
                  P_mf[i_i + j * size_d_index] = P[mf[*i] + j * n];
            }
          else
            //P      = T*(P-K*P(mf,:))*transpose(T)+QQ;
#pragma omp parallel for shared(P_mf)
            for (int i = 0; i < pp; i++)
              for (int j = pure_obs; j < n; j++)
                P_mf[i + j * pp] = P[mf[i] + j * n];

#ifdef BLAS
# pragma omp parallel for shared(K_P)
          for (int i = 0; i < n; i++)
            for (int j = i; j < n; j++)
              {
                double res = 0.0;
                //int j_pp = j * pp;
                for (int k = 0; k < size_d_index; k++)
                  res += K[i + k * n] * P_mf[k + j * size_d_index];
                K_P[i * n + j] = K_P[j * n + i] = res;
              }
          //#pragma omp parallel for shared(P, K_P, P_t_t1)
          for (int i = size_d_index; i < n; i++)
            for (int j = i; j < n; j++)
              {
                unsigned int k = i * n + j;
                P_t_t1[j * n + i] = P_t_t1[k] = P[k] - K_P[k];
              }
          double one = 1.0;
          double zero = 0.0;
          std::copy_n(QQ, n * n, P);
          blas_int n_b = n;
          /*mexPrintf("sizeof(n_b)=%d, n_b=%d, sizeof(n)=%d, n=%d\n",sizeof(n_b),n_b,sizeof(n),n);
            mexEvalString("drawnow;");*/
          dsymm("R", "U", &n_b, &n_b,
                &one, P_t_t1, &n_b,
                T, &n_b, &zero,
                tmp, &n_b);
          dgemm("N", "T", &n_b, &n_b,
                &n_b, &one, tmp, &n_b,
                T, &n_b, &one,
                P, &n_b);
#else
# ifdef CUBLAS
          for (int i = 0; i < n; i++)
            for (int j = i; j < n; j++)
              {
                double res = 0.0;
                //int j_pp = j * pp;
                for (int k = 0; k < size_d_index; k++)
                  res += K[i + k * n] * P_mf[k + j * size_d_index];
                K_P[i * n + j] = K_P[j * n + i] = res;
              }
          //#pragma omp parallel for shared(P, K_P, P_t_t1)
          for (int i = size_d_index; i < n; i++)
            for (int j = i; j < n; j++)
              {
                unsigned int k = i * n + j;
                P_t_t1[j * n + i] = P_t_t1[k] = P[k] - K_P[k];
              }
          mexPrintf("CudaBLAS\n");
          mexEvalString("drawnow;");
          double one = 1.0;
          double zero = 0.0;
          cublasStatus_t status;
          cublasHandle_t handle;
          status = cublasCreate(&handle);
          if (status != CUBLAS_STATUS_SUCCESS)
            {
              mexPrintf("!!!! CUBLAS initialization error\n");
              return false;
            }
          /*int device;
            cudaGetDevice(&device);*/
          int n2 = n * n;
          double *d_A = nullptr, *d_B = nullptr, *d_C = nullptr, *d_D = nullptr;
          // Allocate device memory for the matrices
          if (cudaMalloc(static_cast<void **>(&d_A), n2 * sizeof(double)) != cudaSuccess)
            {
              mexPrintf("!!!! device memory allocation error (allocate A)\n");
              return false;
            }
          if (cudaMalloc(static_cast<void **>(&d_B), n2 * sizeof(d_B[0])) != cudaSuccess)
            {
              mexPrintf("!!!! device memory allocation error (allocate B)\n");
              return false;
            }
          if (cudaMalloc(static_cast<void **>(&d_C), n2 * sizeof(d_C[0])) != cudaSuccess)
            {
              mexPrintf("!!!! device memory allocation error (allocate C)\n");
              return false;
            }
          if (cudaMalloc(static_cast<void **>(&d_D), n2 * sizeof(d_D[0])) != cudaSuccess)
            {
              mexPrintf("!!!! device memory allocation error (allocate D)\n");
              return false;
            }
          // Initialize the device matrices with the host matrices
          status = cublasSetVector(n2, sizeof(P_t_t1[0]), P_t_t1, 1, d_A, 1);
          if (status != CUBLAS_STATUS_SUCCESS)
            {
              mexPrintf("!!!! device access error (write A)\n");
              return false;
            }
          status = cublasSetVector(n2, sizeof(T[0]), T, 1, d_B, 1);
          if (status != CUBLAS_STATUS_SUCCESS)
            {
              mexPrintf("!!!! device access error (write B)\n");
              return false;
            }
          status = cublasSetVector(n2, sizeof(tmp[0]), tmp, 1, d_C, 1);
          if (status != CUBLAS_STATUS_SUCCESS)
            {
              mexPrintf("!!!! device access error (write C)\n");
              return false;
            }
          mexPrintf("just before calling\n");
          mexEvalString("drawnow;");
          status = cublasSetVector(n2, sizeof(QQ[0]), QQ, 1, d_D, 1);
          if (status != CUBLAS_STATUS_SUCCESS)
            {
              mexPrintf("!!!! device access error (write D)\n");
              return false;
            }

          // Performs operation using plain C code

          cublasDsymm(handle, CUBLAS_SIDE_RIGHT, CUBLAS_FILL_MODE_UPPER, n, n,
                      &one, d_A, n,
                      d_B, n, &zero,
                      d_C, n);
          /*dgemm("N", "T", &n_b, &n_b,
            &n_b, &one, tmp, &n_b,
            T, &n_b, &one,
            P, &n_b);*/
          cublasDgemm(handle, CUBLAS_OP_N, CUBLAS_OP_T, n, n,
                      n, &one, d_C, n,
                      d_B, n, &one,
                      d_D, n);
          //double_symm(n, &one, h_A, h_B, &zero, h_C);

          status = cublasGetVector(n2, sizeof(P[0]), d_D, 1, P, 1);
          if (status != CUBLAS_STATUS_SUCCESS)
            {
              mexPrintf("!!!! device access error (read P)\n");
              return false;
            }

# else
#  pragma omp parallel for shared(K_P)
          for (int i = pure_obs; i < n; i++)
            {
              unsigned int i1 = i - pure_obs;
              for (int j = i; j < n; j++)
                {
                  unsigned int j1 = j - pure_obs;
                  double res = 0.0;
                  int j_pp = j * size_d_index;
                  for (int k = 0; k < size_d_index; k++)
                    res += K[i + k * n] * P_mf[k + j_pp];
                  K_P[i1 * n_state + j1] = K_P[j1 * n_state + i1] = res;
                }
            }

#  pragma omp parallel for shared(P_t_t1)
          for (int i = pure_obs; i < n; i++)
            {
              unsigned int i1 = i - pure_obs;
              for (int j = i; j < n; j++)
                {
                  unsigned int j1 = j - pure_obs;
                  unsigned int k1 = i1 * n_state + j1;
                  P_t_t1[j1 * n_state + i1] = P_t_t1[k1] = P[i * n + j] - K_P[k1];
                }
            }

          fill_n(tmp, 0, n * n_state);

#  pragma omp parallel for shared(tmp)
          for (int i = 0; i < n; i++)
            {
              int max_k = i_nz_state_var[i];
              for (int j = pure_obs; j < n; j++)
                {
                  int j1 = j - pure_obs;
                  int j1_n_state = j1 * n_state - pure_obs;
                  int indx_tmp = i + j1 * n;
                  for (int k = pure_obs; k < max_k; k++)
                    tmp[indx_tmp] += T[i + k * n] * P_t_t1[k + j1_n_state];
                }
            }

          fill_n(P, 0, n * n);

          int n_n_obs = -n * pure_obs;
#  pragma omp parallel for shared(P)
          for (int i = 0; i < n; i++)
            {
              for (int j = i; j < n; j++)
                {
                  int max_k = i_nz_state_var[j];
                  int P_indx = i * n + j;
                  for (int k = pure_obs; k < max_k; k++)
                    {
                      int k_n = k * n;
                      P[P_indx] += tmp[i + k_n + n_n_obs] * T[j + k_n];
                    }
                }
            }

#  pragma omp parallel for shared(P)
          for (int i = 0; i < n; i++)
            {
              for (int j = i; j < n; j++)
                P[j + i * n] += QQ[j + i * n];
              for (int j = i + 1; j < n; j++)
                P[i + j * n] = P[j + i * n];
            }
# endif
#endif
          if (t >= no_more_missing_observations)
            {
              double max_abs = 0.0;
              for (int i = 0; i < n * size_d_index; i++)
                {
                  double res = abs(K[i] - oldK[i]);
                  max_abs = std::max(res, max_abs);
                }
              notsteady = max_abs > riccati_tol;

              //oldK = K(:);

              std::copy_n(K, n * pp, oldK.get());
            }
        }
      t++;
    }

  if (F_singular)
    mexErrMsgTxt("The variance of the forecast error remains singular until the end of the sample\n");
  if (t < smpl)
    block_kalman_filter_ss();
  return true;
}

void
BlockKalmanFilter::return_results_and_clean(int nlhs, mxArray *plhs[])
{
  if (nlhs > 2)
    mexErrMsgTxt("block_kalman_filter provides at most 2 output argument.");

  if (nlhs >= 1)
    {
      plhs[0] = mxCreateDoubleMatrix(1, 1, mxREAL);
      double *pind = mxGetPr(plhs[0]);
      pind[0] = LIK;
    }

  if (nlhs == 2)
    plhs[1] = plik;
  else
    mxDestroyArray(plik);

  mxDestroyArray(pa);
  mxDestroyArray(p_tmp);
  mxDestroyArray(pQQ);
  mxDestroyArray(pv);
  mxDestroyArray(pF);
  mxDestroyArray(piF);
  mxDestroyArray(p_P_t_t1);
  mxDestroyArray(pK);
  mxDestroyArray(p_K_P);
}

void
mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
  BlockKalmanFilter block_kalman_filter(nrhs, prhs);
  if (block_kalman_filter.block_kalman_filter(nlhs, plhs))
    block_kalman_filter.return_results_and_clean(nlhs, plhs);
}
