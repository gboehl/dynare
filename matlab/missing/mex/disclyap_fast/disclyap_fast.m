function [X,exitflag]=disclyap_fast(G,V,tol,check_flag,max_iter)
% [X,exitflag]=disclyap_fast(G,V,tol,check_flag)
% Inputs:
%   - G             [double]    first input matrix
%   - V             [double]    second input matrix
%   - tol           [scalar]    tolerance criterion
%   - check_flag    [boolean]   if true: check positive-definiteness
%   - max_iter      [scalar]    maximum number of iterations

% Outputs:
%   - X             [double]    solution matrix
%   - exitflag      [scalar]    0 if solution is found, 1 otherwise
%
% Solve the discrete Lyapunov Equation
% X=G*X*G'+V
% Using the Doubling Algorithm
%
% If check_flag is defined then the code will check if the resulting X
% is positive definite and generate an error message if it is not
%
% Joe Pearlman and Alejandro Justiniano
% 3/5/2005

% Copyright (C) 2010-2020 Dynare Team
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

if nargin <= 3 || isempty(check_flag)
    check_flag = 0;
end
if nargin <=4
    max_iter=2000;
end
exitflag=0;

P0=V;
A0=G;

matd=1;
iter=1;
while matd > tol && iter< max_iter
    P1=P0+A0*P0*A0';
    A1=A0*A0;
    matd=max( max( abs( P1 - P0 ) ) );
    P0=P1;
    A0=A1;
    iter=iter+1;
end
if iter==max_iter || isnan(matd)
    X=NaN(size(P0));
    exitflag=1;
    return
end

X=(P0+P0')/2;

% Check that X is positive definite
if check_flag==1
    [~,p]=chol(X);
    if p ~= 0
        exitflag=1;
    end
end