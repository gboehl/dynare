function [ss, tt, zz, sdim, eigval, info] = mjdgges(e, d, qz_criterium, zhreshold) % --*-- Unitary tests --*--
%
% INPUTS
%   e            [double] real square (n*n) matrix.
%   d            [double] real square (n*n) matrix.
%   qz_criterium [double] scalar (1+epsilon).
%   zhreshold    [double] ignored (in the DLL, this parameter is used for
%                                  detecting eigenvalues too close to 0/0)
%
% OUTPUTS
%   ss           [double] (n*n) quasi-triangular matrix.
%   tt           [double] (n*n) quasi-triangular matrix.
%   zz           [double] (n*n) orthogonal matrix.
%   sdim         [integer] scalar: number of stable eigenvalues.
%   eigval       [complex] (n*1) vector of generalized eigenvalues.
%   info         [integer] scalar.
%
% SPECIAL REQUIREMENTS
%   none.

% Copyright (C) 1996-2020 Dynare Team
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

if nargin > 5 || nargin < 2 || nargout > 7 || nargout == 0
    error('MJDGGES: takes 2, 3 or 4 input arguments and between 1 and 7 output arguments.')
end

if isoctave
    error('Octave unsupported, since it does not have real qz, ordqz and ordeig')
end

[me, ne] = size(e);
[md, nd] = size(d);
if ~isreal(e) || ~isreal(d) || me ~= ne || md ~= nd || me ~= nd
    error('MJDGGES requires two square real matrices of the same dimension.')
end

if nargin < 3
    qz_criterium = 1 + 1e-6;
end

info = 0;

try
    [ss, tt, qq, zz] = qz(e, d, 'real');
    eigval = ordeig(ss, tt);
    select = abs(eigval) < qz_criterium;
    sdim = sum(select);
    [ss, tt, qq, zz] = ordqz(ss, tt, qq, zz, select);
    eigval = ordeig(ss, tt);
catch
    info = 1; % Not as precise as lapack's info!
end

%@test:1
%$ try
%$     E =[0,0,0,0,0,0,0,-1,0,0;0,0,0,0,0,0,0,0,-1,0;0,0,0,0,0,0,0,0,0,-0.990099009900990;0,0,0,0,0,0,0,0,0,0.0990099009900990;0,0,0,-1.01010101010101,0,0.0427672955974843,0,0,0,0;0,0,0,0,0,0.128301886792453,-1,0,0,0;0.800000000000000,0,0,0,0,0,0,0,0,0;0,1,0,0,0,0,1,0,0,0;0,0,0.900000000000000,0,0,0,0,0,0,0;0,0,0,-1.01010101010101,-1,0,2,0,-1,0];
%$     D=[0,0,0,0,-1,0,0,-0.792000000000000,0,0;0,0,0,0,0,0,0,0,-0.990000000000000,0;0,0,0,-0.000493818030899887,0,0,0,0,0,-0.882178217821782;0,0,0,-1.00493818030900,0,0,0,0,0,0.0882178217821782;0,0,0,-1,0.128301886792453,0,0,0,0,0;-1,0,0,0,0,0,-0.990000000000000,0,0,0;1,0,0,0,0,0,0,0,0,0;0,1,0,0,0,0,0,0,0,0;0,0,1,0,0,0,0,0,0,0;0,0,0,0,-1,0,0,0,0,0];
%$     [ss, tt, w, sdim, dr.eigval, info1]=mjdgges(E, D, 1.000001, 1e-06);
%$     if sdim==5
%$         t(1) = 1;
%$     else 
%$         t(1) = 0;
%$     end
%$ catch
%$     t(1) = 0;
%$ end
%$ T = all(t);
%@eof:1