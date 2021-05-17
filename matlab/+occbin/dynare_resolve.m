function [A,B,ys,info,M_,oo_,TT, RR, CC, A0, B0] ...
    = dynare_resolve(M_,options_,oo_,regime_history, reduced_state_space, A, B)
% function [A,B,ys,info,M_,options_,oo_,TT, RR, CC, A0, B0] ...
%     = dynare_resolve(M_,options_,oo_,regime_history, reduced_state_space, A, B)
% Computes the linear approximation and the matrices A and B of the
% transition equation. Mirrors dynare_resolve
%
% Inputs:
% - M_                  [structure]     Matlab's structure describing the model
% - options_            [structure]     Matlab's structure containing the options
% - oo_                 [structure]     Matlab's structure containing the results
% - reduced_state_space [string]
% - A                   [double]        State transition matrix
% - B                   [double]        shock impact matrix
% Outputs:
% - A                   [double]        State transition matrix (potentially for restricted state space)
% - B                   [double]        shock impact matrix (potentially for restricted state space)
% - ys                  [double]        vector of steady state values
% - info                [double]        4 by 1 vector with exit flag and information
% - M_                  [structure]     Matlab's structure describing the model
% - options_            [structure]     Matlab's structure containing the options
% - oo_                 [structure]     Matlab's structure containing the results
% - TT                  [N by N]            transition matrix of state space for each period
% - RR                  [N by N_exo by T]   shock impact matrix of state space for each period
% - CC                  [N by N_exo by T]   constant of state space for each period
% - A0                  [double]        State transition matrix (unrestricted state space)
% - B0                  [double]        shock impact matrix (unrestricted state space)

% Copyright (C) 2001-2021 Dynare Team
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

if nargin<6
    [A,B,ys,info,M_,oo_] = dynare_resolve(M_,options_,oo_);
else
    ys = oo_.dr.ys;
    info = 0;
end
if  ~info(1) && nargin>4 && ~isempty(regime_history)
    opts_regime.regime_history=regime_history;
    opts_regime.binding_indicator=[];
    [TT, RR, CC] = ...
        occbin.check_regimes(A, B, [], opts_regime, M_,oo_,options_);
else
    TT=A;
    RR=B;
    CC=zeros(size(A,1),1);
end

A0=A;
B0=B;
if  ~info(1) && nargin>4 && ischar(reduced_state_space) && ~isempty(reduced_state_space)
    iv = oo_.dr.restrict_var_list;
    A=A(iv,iv);
    B=B(iv,:);
    TT=TT(iv,iv,:);
    RR=RR(iv,:,:);
    CC=CC(iv,:);
end