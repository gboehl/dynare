/*
 * Copyright © 2004-2011 Ondra Kamenik
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

#include "SylvParams.hh"

#include <array>
#include <iostream>

void
SylvParams::print(const std::string& prefix) const
{
  print(std::cout, prefix);
}

void
SylvParams::print(std::ostream& fdesc, const std::string& prefix) const
{
  method.print(fdesc, prefix, "method             ");
  rcondA1.print(fdesc, prefix, "reci. cond1 A      ");
  rcondAI.print(fdesc, prefix, "reci. cond∞ A      ");
  bs_norm.print(fdesc, prefix, "log₁₀ diag norm    ");
  f_err1.print(fdesc, prefix, "abs. err 1 F diag  ");
  f_errI.print(fdesc, prefix, "abs. err ∞ F diag  ");
  viv_err1.print(fdesc, prefix, "abs. err 1 V·V⁻¹   ");
  viv_errI.print(fdesc, prefix, "abs. err ∞ V·V⁻¹   ");
  ivv_err1.print(fdesc, prefix, "abs. err 1 V⁻¹·V   ");
  ivv_errI.print(fdesc, prefix, "abs. err ∞ V⁻¹·V   ");
  f_blocks.print(fdesc, prefix, "num blocks in F    ");
  f_largest.print(fdesc, prefix, "largest block in F ");
  f_zeros.print(fdesc, prefix, "num zeros in F     ");
  f_offdiag.print(fdesc, prefix, "num offdiag in F   ");
  if (*method == solve_method::iter)
    {
      converged.print(fdesc, prefix, "converged          ");
      convergence_tol.print(fdesc, prefix, "convergence tol.   ");
      iter_last_norm.print(fdesc, prefix, "last norm          ");
      max_num_iter.print(fdesc, prefix, "max num iter       ");
      num_iter.print(fdesc, prefix, "num iter           ");
    }
  else
    eig_min.print(fdesc, prefix, "minimum eigenvalue ");

  mat_err1.print(fdesc, prefix, "rel. matrix norm1  ");
  mat_errI.print(fdesc, prefix, "rel. matrix norm∞  ");
  mat_errF.print(fdesc, prefix, "rel. matrix normFro");
  vec_err1.print(fdesc, prefix, "rel. vector norm1  ");
  vec_errI.print(fdesc, prefix, "rel. vector norm∞  ");
  cpu_time.print(fdesc, prefix, "time (CPU secs)    ");
}

std::vector<const char*>
SylvParams::getArrayNames() const
{
  std::vector<const char*> names;
  if (method.getStatus() != status::undef)
    names.push_back("method");
  if (convergence_tol.getStatus() != status::undef)
    names.push_back("convergence_tol");
  if (max_num_iter.getStatus() != status::undef)
    names.push_back("max_num_iter");
  if (bs_norm.getStatus() != status::undef)
    names.push_back("bs_norm");
  if (converged.getStatus() != status::undef)
    names.push_back("converged");
  if (iter_last_norm.getStatus() != status::undef)
    names.push_back("iter_last_norm");
  if (num_iter.getStatus() != status::undef)
    names.push_back("num_iter");
  if (f_err1.getStatus() != status::undef)
    names.push_back("f_err1");
  if (f_errI.getStatus() != status::undef)
    names.push_back("f_errI");
  if (viv_err1.getStatus() != status::undef)
    names.push_back("viv_err1");
  if (viv_errI.getStatus() != status::undef)
    names.push_back("viv_errI");
  if (ivv_err1.getStatus() != status::undef)
    names.push_back("ivv_err1");
  if (ivv_errI.getStatus() != status::undef)
    names.push_back("ivv_errI");
  if (f_blocks.getStatus() != status::undef)
    names.push_back("f_blocks");
  if (f_largest.getStatus() != status::undef)
    names.push_back("f_largest");
  if (f_zeros.getStatus() != status::undef)
    names.push_back("f_zeros");
  if (f_offdiag.getStatus() != status::undef)
    names.push_back("f_offdiag");
  if (rcondA1.getStatus() != status::undef)
    names.push_back("rcondA1");
  if (rcondAI.getStatus() != status::undef)
    names.push_back("rcondAI");
  if (eig_min.getStatus() != status::undef)
    names.push_back("eig_min");
  if (mat_err1.getStatus() != status::undef)
    names.push_back("mat_err1");
  if (mat_errI.getStatus() != status::undef)
    names.push_back("mat_errI");
  if (mat_errF.getStatus() != status::undef)
    names.push_back("mat_errF");
  if (vec_err1.getStatus() != status::undef)
    names.push_back("vec_err1");
  if (vec_errI.getStatus() != status::undef)
    names.push_back("vec_errI");
  if (cpu_time.getStatus() != status::undef)
    names.push_back("cpu_time");

  return names;
}

mxArray*
SylvParams::DoubleParamItem::createMatlabArray() const
{
  return mxCreateDoubleScalar(value);
}

mxArray*
SylvParams::IntParamItem::createMatlabArray() const
{
  mxArray* res = mxCreateNumericMatrix(1, 1, mxINT32_CLASS, mxREAL);
#if MX_HAS_INTERLEAVED_COMPLEX
  *mxGetInt32s(res) = value;
#else
  *static_cast<int*>(mxGetData(res)) = value;
#endif
  return res;
}

mxArray*
SylvParams::BoolParamItem::createMatlabArray() const
{
  if (value)
    return mxCreateString("true");
  else
    return mxCreateString("false");
}

mxArray*
SylvParams::MethodParamItem::createMatlabArray() const
{
  if (value == solve_method::iter)
    return mxCreateString("iterative");
  else
    return mxCreateString("recursive");
}

mxArray*
SylvParams::createStructArray() const
{
  auto names = getArrayNames();
  const std::array<mwSize, 2> dims {1, 1};
  mxArray* const res = mxCreateStructArray(dims.size(), dims.data(), names.size(), names.data());

  int i = 0;
  if (method.getStatus() != status::undef)
    mxSetFieldByNumber(res, 0, i++, method.createMatlabArray());
  if (convergence_tol.getStatus() != status::undef)
    mxSetFieldByNumber(res, 0, i++, convergence_tol.createMatlabArray());
  if (max_num_iter.getStatus() != status::undef)
    mxSetFieldByNumber(res, 0, i++, max_num_iter.createMatlabArray());
  if (bs_norm.getStatus() != status::undef)
    mxSetFieldByNumber(res, 0, i++, bs_norm.createMatlabArray());
  if (converged.getStatus() != status::undef)
    mxSetFieldByNumber(res, 0, i++, converged.createMatlabArray());
  if (iter_last_norm.getStatus() != status::undef)
    mxSetFieldByNumber(res, 0, i++, iter_last_norm.createMatlabArray());
  if (num_iter.getStatus() != status::undef)
    mxSetFieldByNumber(res, 0, i++, num_iter.createMatlabArray());
  if (f_err1.getStatus() != status::undef)
    mxSetFieldByNumber(res, 0, i++, f_err1.createMatlabArray());
  if (f_errI.getStatus() != status::undef)
    mxSetFieldByNumber(res, 0, i++, f_errI.createMatlabArray());
  if (viv_err1.getStatus() != status::undef)
    mxSetFieldByNumber(res, 0, i++, viv_err1.createMatlabArray());
  if (viv_errI.getStatus() != status::undef)
    mxSetFieldByNumber(res, 0, i++, viv_errI.createMatlabArray());
  if (ivv_err1.getStatus() != status::undef)
    mxSetFieldByNumber(res, 0, i++, ivv_err1.createMatlabArray());
  if (ivv_errI.getStatus() != status::undef)
    mxSetFieldByNumber(res, 0, i++, ivv_errI.createMatlabArray());
  if (f_blocks.getStatus() != status::undef)
    mxSetFieldByNumber(res, 0, i++, f_blocks.createMatlabArray());
  if (f_largest.getStatus() != status::undef)
    mxSetFieldByNumber(res, 0, i++, f_largest.createMatlabArray());
  if (f_zeros.getStatus() != status::undef)
    mxSetFieldByNumber(res, 0, i++, f_zeros.createMatlabArray());
  if (f_offdiag.getStatus() != status::undef)
    mxSetFieldByNumber(res, 0, i++, f_offdiag.createMatlabArray());
  if (rcondA1.getStatus() != status::undef)
    mxSetFieldByNumber(res, 0, i++, rcondA1.createMatlabArray());
  if (rcondAI.getStatus() != status::undef)
    mxSetFieldByNumber(res, 0, i++, rcondAI.createMatlabArray());
  if (eig_min.getStatus() != status::undef)
    mxSetFieldByNumber(res, 0, i++, eig_min.createMatlabArray());
  if (mat_err1.getStatus() != status::undef)
    mxSetFieldByNumber(res, 0, i++, mat_err1.createMatlabArray());
  if (mat_errI.getStatus() != status::undef)
    mxSetFieldByNumber(res, 0, i++, mat_errI.createMatlabArray());
  if (mat_errF.getStatus() != status::undef)
    mxSetFieldByNumber(res, 0, i++, mat_errF.createMatlabArray());
  if (vec_err1.getStatus() != status::undef)
    mxSetFieldByNumber(res, 0, i++, vec_err1.createMatlabArray());
  if (vec_errI.getStatus() != status::undef)
    mxSetFieldByNumber(res, 0, i++, vec_errI.createMatlabArray());
  if (cpu_time.getStatus() != status::undef)
    mxSetFieldByNumber(res, 0, i++, cpu_time.createMatlabArray());

  return res;
}
