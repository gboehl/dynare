function [h0, h1, longruncorrectionparameter] = hVectors(params, H, ids, idns, auxmodel)

% INPUTS
% - params          [double]     (m+1)*1 vector, PAC parameters (lag polynomial coefficients and discount factor).
% - H               [double]     (np*np) matrix, companion matrix of the VAR(p) model.
% - ids             [integer]    n*1 vector, selection of the stationary components of the VAR.
% - idns            [integer]    n*1 vector, selection of the non stationary components of the VAR.
%
% OUTPUTS
% - h0              [double]     1*n vector.
% - h1              [double]     1*n vector.
%
% REMARKS
% The parameters are ordered as follows in the params vector:
%
%    params(1)       ⟶ Error correction parameter.
%    params(2:end-1) ⟶ Autoregressive parameters.
%    params(end)     ⟶ Discount factor.

% Copyright (C) 2018 Dynare Team
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

if nargout>1 && nargin<4
    error('Wrong number of Inputs/Outputs!')
end

[G, alpha, beta] = buildGmatrixWithAlphaAndBeta(params);

A = [alpha; 1];

A_1 = polyval(A, 1.0);
A_b = polyval(A, beta);

m = length(alpha);
n = length(H);

tmp = eye(n*m)-kron(G, transpose(H));

if isempty(ids)
    h0 = [];
else
    h0 = A_1*A_b*((kron(iota(m, m), H))'*(tmp\kron(iota(m, m), iota(n, ids))));
end

if nargout>1
    if isempty(idns)
        h1 = [];
    else
        switch auxmodel
          case {'var', 'trend_component'}
            h1 = A_1*A_b*(kron(iota(m, m)'*inv(eye(m)-G), H')*(tmp\kron(iota(m, m), iota(n, idns))));
          case 'Need to check in which case we should trigger this one...'
            h1 = A_1*A_b*(kron(iota(m, m)'*inv(eye(m)-G), (H'-eye(length(H))))*(tmp\kron(iota(m, m), iota(n, idns))));
        end
    end
end

if nargout>2
    d = A_1*A_b*(iota(m, m)'*inv((eye(m)-G)*(eye(m)-G))*iota(m, m));
    longruncorrectionparameter = (1-sum(params(2:end-1))-d);
end