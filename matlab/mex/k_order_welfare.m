% [unc_welfare, U_dynpp_derivs, W_dynpp_derivs, U_dyn_derivs, W_dyn_derivs] = k_order_welfare(dr, DynareModel, DynareOptions)
% computes a k-th order approximation of welfare
%
% INPUTS
% dr:              struct   describing the reduced form solution of the model.
% DynareModel:     struct   jobs's parameters
% DynareOptions:   struct   job's options
%
% OUTPUTS
%
% unc_welfare      double   approximation of conditional welfare from the non-stochastic steady state and allowing stochastic shocks
%
% U_dynpp_derivs   struct   Derivatives of the felicity function in Dynare++ format.
%                         The tensors are folded and the Taylor coefficients (1/n!) are
%                         included.
%                         Fieldnames are U_0, U_1, U_2, …
% W_dynpp_derivs   struct   Derivatives of the welfare function in Dynare++ format.
%                         The tensors are folded and the Taylor coefficients (1/n!) are
%                         included.
%                         Fieldnames are W_0, W_1, W_2, …
% U_dyn_derivs     struct   Derivatives of the felicity function in Dynare format, up to third order.
%                         Matrices are dense and unfolded. The Taylor coefficients (1/2
%                         and 1/6) aren't included.
%                         The derivatives w.r.t. different categories of variables
%                         are separated.
%                         Fieldnames are:
%                          + Uy, Uu
%                          + if order ≥ 2: Uyy, Uyu, Uuu, Uss
%                          + if order ≥ 3: Uyyy, Uyyu, Uyuu, Uuuu, Uyss, Uuss
%
% W_dyn_derivs     struct   Derivatives of the welfare function in Dynare format, up to third order.
%                         Matrices are dense and unfolded. The Taylor coefficients (1/2
%                         and 1/6) aren't included.
%                         The derivatives w.r.t. different categories of variables
%                         are separated.
%                         Fieldnames are:
%                          + Wy, Wu
%                          + if order ≥ 2: Wyy, Wyu, Wuu, Wss
%                          + if order ≥ 3: Wyyy, Wyyu, Wyuu, Wuuu, Wyss, Wuss
%
% k_order_welfare is a compiled MEX function. Its source code is in
% dynare/mex/sources/k_order_welfare/k_order_welfare.cc and it uses code provided by
% dynare++

% Copyright (C) 2021 Dynare Team
%
% This file is part of Dynare.
%
% Dynare is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
%
% Dynare is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with Dynare.  If not, see <https://www.gnu.org/licenses/>.
