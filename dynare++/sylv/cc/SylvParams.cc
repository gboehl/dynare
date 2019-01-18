/* $Header: /var/lib/cvs/dynare_cpp/sylv/cc/SylvParams.cpp,v 1.1.1.1 2004/06/04 13:00:52 kamenik Exp $ */

/* Tag $Name:  $ */

#include "SylvParams.hh"

#include <iostream>

void
SylvParams::print(const std::string &prefix) const
{
  print(std::cout, prefix);
}

void
SylvParams::print(std::ostream &fdesc, const std::string &prefix) const
{
  method.print(fdesc, prefix, "method             ");
  rcondA1.print(fdesc, prefix,  "reci. cond1 A      ");
  rcondAI.print(fdesc, prefix,  "reci. condInf A    ");
  bs_norm.print(fdesc, prefix,  "log10 diag norm    ");
  f_err1.print(fdesc, prefix,   "abs. err 1 F diag  ");
  f_errI.print(fdesc, prefix,   "abs. err I F diag  ");
  viv_err1.print(fdesc, prefix, "abs. err 1 V*invV  ");
  viv_errI.print(fdesc, prefix, "abs. err I V*invV  ");
  ivv_err1.print(fdesc, prefix, "abs. err 1 invV*V  ");
  ivv_errI.print(fdesc, prefix, "abs. err I invV*V  ");
  f_blocks.print(fdesc, prefix, "num blocks in F    ");
  f_largest.print(fdesc, prefix, "largest block in F ");
  f_zeros.print(fdesc, prefix,  "num zeros in F     ");
  f_offdiag.print(fdesc, prefix, "num offdiag in F   ");
  if (*method == iter)
    {
      converged.print(fdesc, prefix,       "converged          ");
      convergence_tol.print(fdesc, prefix, "convergence tol.   ");
      iter_last_norm.print(fdesc, prefix,  "last norm          ");
      max_num_iter.print(fdesc, prefix,    "max num iter       ");
      num_iter.print(fdesc, prefix,        "num iter           ");
    }
  else
    eig_min.print(fdesc, prefix,         "minimum eigenvalue ");

  mat_err1.print(fdesc, prefix, "rel. matrix norm1  ");
  mat_errI.print(fdesc, prefix, "rel. matrix normInf");
  mat_errF.print(fdesc, prefix, "rel. matrix normFro");
  vec_err1.print(fdesc, prefix, "rel. vector norm1  ");
  vec_errI.print(fdesc, prefix, "rel. vector normInf");
  cpu_time.print(fdesc, prefix, "time (CPU secs)    ");
}

void
SylvParams::setArrayNames(int &num, const char **names) const
{
  num = 0;
  if (method.getStatus() != undef)
    names[num++] = "method";
  if (convergence_tol.getStatus() != undef)
    names[num++] = "convergence_tol";
  if (max_num_iter.getStatus() != undef)
    names[num++] = "max_num_iter";
  if (bs_norm.getStatus() != undef)
    names[num++] = "bs_norm";
  if (converged.getStatus() != undef)
    names[num++] = "converged";
  if (iter_last_norm.getStatus() != undef)
    names[num++] = "iter_last_norm";
  if (num_iter.getStatus() != undef)
    names[num++] = "num_iter";
  if (f_err1.getStatus() != undef)
    names[num++] = "f_err1";
  if (f_errI.getStatus() != undef)
    names[num++] = "f_errI";
  if (viv_err1.getStatus() != undef)
    names[num++] = "viv_err1";
  if (viv_errI.getStatus() != undef)
    names[num++] = "viv_errI";
  if (ivv_err1.getStatus() != undef)
    names[num++] = "ivv_err1";
  if (ivv_errI.getStatus() != undef)
    names[num++] = "ivv_errI";
  if (f_blocks.getStatus() != undef)
    names[num++] = "f_blocks";
  if (f_largest.getStatus() != undef)
    names[num++] = "f_largest";
  if (f_zeros.getStatus() != undef)
    names[num++] = "f_zeros";
  if (f_offdiag.getStatus() != undef)
    names[num++] = "f_offdiag";
  if (rcondA1.getStatus() != undef)
    names[num++] = "rcondA1";
  if (rcondAI.getStatus() != undef)
    names[num++] = "rcondAI";
  if (eig_min.getStatus() != undef)
    names[num++] = "eig_min";
  if (mat_err1.getStatus() != undef)
    names[num++] = "mat_err1";
  if (mat_errI.getStatus() != undef)
    names[num++] = "mat_errI";
  if (mat_errF.getStatus() != undef)
    names[num++] = "mat_errF";
  if (vec_err1.getStatus() != undef)
    names[num++] = "vec_err1";
  if (vec_errI.getStatus() != undef)
    names[num++] = "vec_errI";
  if (cpu_time.getStatus() != undef)
    names[num++] = "cpu_time";
}

#if defined(MATLAB_MEX_FILE) || defined(OCTAVE_MEX_FILE)
mxArray *
SylvParams::DoubleParamItem::createMatlabArray() const
{
  return mxCreateDoubleScalar(value);
}

mxArray *
SylvParams::IntParamItem::createMatlabArray() const
{
  mxArray *res = mxCreateNumericMatrix(1, 1, mxINT32_CLASS, mxREAL);
  *((int *) mxGetData(res)) = value;
  return res;
}

mxArray *
SylvParams::BoolParamItem::createMatlabArray() const
{
  if (value)
    return mxCreateString("true");
  else
    return mxCreateString("false");
}

mxArray *
SylvParams::MethodParamItem::createMatlabArray() const
{
  if (value == iter)
    return mxCreateString("iterative");
  else
    return mxCreateString("recursive");
}

mxArray *
SylvParams::createStructArray() const
{
  const char *names[50];
  int num;
  setArrayNames(num, names);
  const mwSize dims[] = {1, 1};
  mxArray *const res = mxCreateStructArray(2, dims, num, names);

  int i = 0;
  if (method.getStatus() != undef)
    mxSetFieldByNumber(res, 0, i++, method.createMatlabArray());
  if (convergence_tol.getStatus() != undef)
    mxSetFieldByNumber(res, 0, i++, convergence_tol.createMatlabArray());
  if (max_num_iter.getStatus() != undef)
    mxSetFieldByNumber(res, 0, i++, max_num_iter.createMatlabArray());
  if (bs_norm.getStatus() != undef)
    mxSetFieldByNumber(res, 0, i++, bs_norm.createMatlabArray());
  if (converged.getStatus() != undef)
    mxSetFieldByNumber(res, 0, i++, converged.createMatlabArray());
  if (iter_last_norm.getStatus() != undef)
    mxSetFieldByNumber(res, 0, i++, iter_last_norm.createMatlabArray());
  if (num_iter.getStatus() != undef)
    mxSetFieldByNumber(res, 0, i++, num_iter.createMatlabArray());
  if (f_err1.getStatus() != undef)
    mxSetFieldByNumber(res, 0, i++, f_err1.createMatlabArray());
  if (f_errI.getStatus() != undef)
    mxSetFieldByNumber(res, 0, i++, f_errI.createMatlabArray());
  if (viv_err1.getStatus() != undef)
    mxSetFieldByNumber(res, 0, i++, viv_err1.createMatlabArray());
  if (viv_errI.getStatus() != undef)
    mxSetFieldByNumber(res, 0, i++, viv_errI.createMatlabArray());
  if (ivv_err1.getStatus() != undef)
    mxSetFieldByNumber(res, 0, i++, ivv_err1.createMatlabArray());
  if (ivv_errI.getStatus() != undef)
    mxSetFieldByNumber(res, 0, i++, ivv_errI.createMatlabArray());
  if (f_blocks.getStatus() != undef)
    mxSetFieldByNumber(res, 0, i++, f_blocks.createMatlabArray());
  if (f_largest.getStatus() != undef)
    mxSetFieldByNumber(res, 0, i++, f_largest.createMatlabArray());
  if (f_zeros.getStatus() != undef)
    mxSetFieldByNumber(res, 0, i++, f_zeros.createMatlabArray());
  if (f_offdiag.getStatus() != undef)
    mxSetFieldByNumber(res, 0, i++, f_offdiag.createMatlabArray());
  if (rcondA1.getStatus() != undef)
    mxSetFieldByNumber(res, 0, i++, rcondA1.createMatlabArray());
  if (rcondAI.getStatus() != undef)
    mxSetFieldByNumber(res, 0, i++, rcondAI.createMatlabArray());
  if (eig_min.getStatus() != undef)
    mxSetFieldByNumber(res, 0, i++, eig_min.createMatlabArray());
  if (mat_err1.getStatus() != undef)
    mxSetFieldByNumber(res, 0, i++, mat_err1.createMatlabArray());
  if (mat_errI.getStatus() != undef)
    mxSetFieldByNumber(res, 0, i++, mat_errI.createMatlabArray());
  if (mat_errF.getStatus() != undef)
    mxSetFieldByNumber(res, 0, i++, mat_errF.createMatlabArray());
  if (vec_err1.getStatus() != undef)
    mxSetFieldByNumber(res, 0, i++, vec_err1.createMatlabArray());
  if (vec_errI.getStatus() != undef)
    mxSetFieldByNumber(res, 0, i++, vec_errI.createMatlabArray());
  if (cpu_time.getStatus() != undef)
    mxSetFieldByNumber(res, 0, i++, cpu_time.createMatlabArray());

  return res;
}
#endif
