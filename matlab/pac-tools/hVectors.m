function [h, lrcp] = hVectors(params, H, auxmodel, kind, id)

% INPUTS
% - params          [double]     (m+1)*1 vector, PAC parameters (lag polynomial coefficients and discount factor).
% - H               [double]     (np*np) matrix, companion matrix of the VAR(p) model.
% - auxmodel        [char]       kind of auxiliary model, possible values are: 'var' and 'trend_component'.
% - kind            [char]       kind of expectation in PAC equation (See FRB/US doc), possible values are 'll', 'dl' and 'dd'.
% - id              [integer]    scalar, index pointing to the target in the auxiliary model.
%
% OUTPUTS
% - h               [double]     1*n vector of weights (used to compute the linear combination of the variables in the companion VAR representation of the auxiliary model).
% - lrcp            [double]     scalar, parameter for the growth neutrality correction.
%
% REMARKS
% The parameters are ordered as follows in the params vector:
%
%    params(1)       ⟶ Error correction parameter.
%    params(2:end-1) ⟶ Autoregressive parameters.
%    params(end)     ⟶ Discount factor.

% Copyright © 2018-2021 Dynare Team
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

tmp = eye(n*m)-kron(G, transpose(H)); % inv(W2)

switch kind
  case 'll'
    h = A_1*A_b*((kron(iota(m, m), H))'*(tmp\kron(iota(m, m), iota(n, id))));
  case 'dd'
    h = A_1*A_b*(kron(iota(m, m)'*inv(eye(m)-G), H')*(tmp\kron(iota(m, m), iota(n, id))));
  case 'dl'
    h = A_1*A_b*(kron(iota(m, m)'*inv(eye(m)-G), (H'-eye(length(H))))*(tmp\kron(iota(m, m), iota(n, id))));
  otherwise
    error('Unknown kind value in PAC model.')
end

if nargin>1
    if isequal(kind, 'll')
        lrcp = NaN;
    else
        d = A_1*A_b*(iota(m, m)'*inv((eye(m)-G)*(eye(m)-G))*iota(m, m));
        lrcp = (1-sum(params(2:end-1))-d);
    end
end