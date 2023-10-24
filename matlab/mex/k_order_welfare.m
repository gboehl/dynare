% W = k_order_welfare(dr, M_, options_)
% computes a k-th order approximation of welfare
%
% INPUTS
% dr:              struct   describing the reduced form solution of the model.
% M_:              struct   jobs's parameters
% options_:        struct   job's options
%
% OUTPUTS
%
% W                struct   Derivatives of the welfare function in Dynare++ format.
%                         The tensors are folded and the Taylor coefficients (1/n!) are
%                         included.
%                         Fieldnames are W_0, W_1, W_2, …
%
% k_order_welfare is a compiled MEX function. Its source code is in
% dynare/mex/sources/k_order_welfare/k_order_welfare.cc and it uses code provided by
% dynare++

% Copyright © 2021-2022 Dynare Team
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
