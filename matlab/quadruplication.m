function [QP,QPinv] = quadruplication(p)
% Computes the Quadruplication Matrix QP (and its Moore-Penrose inverse)
% such that for any p-dimensional vector x:
% y=kron(kron(kron(x,x),x),x)=QP*z 
% where z is of dimension np=p*(p+1)*(p+2)*(p+3)/2 and is obtained from y 
% by removing each second and later occurence of the same element.
% This is a generalization of the Duplication matrix.
% Reference: Meijer (2005) - Matrix algebra for higher order moments.
%            Linear Algebra and its Applications, 410,pp. 112-134
% =========================================================================
% INPUTS
%    * p        [integer]    size of vector
% -------------------------------------------------------------------------
% OUTPUTS
%    * QP       [p^4 by np]  Quadruplication matrix
%    * QPinv    [np by np]   Moore-Penrose inverse of QP
% -------------------------------------------------------------------------
% This function is called by
%   * pruned_state_space_system.m
% -------------------------------------------------------------------------
% This function calls
%   * mue (embedded)
%   * uperm
% =========================================================================
% Copyright © 2020 Dynare Team
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
% =========================================================================
np = p*(p+1)*(p+2)*(p+3)/24;
QP = spalloc(p^4,np,p^4);
if nargout > 1
    QPinv = spalloc(np,np,p^4);
end
counti = 1;
for l=1:p
    for k=l:p
        for j=k:p
            for i=j:p                             
                idx = uperm([i j k l]);
                for r = 1:size(idx,1)
                    ii = idx(r,1); jj= idx(r,2); kk=idx(r,3); ll=idx(r,4);
                    n = ii + (jj-1)*p + (kk-1)*p^2 + (ll-1)*p^3;                    
                    m = mue(p,i,j,k,l);
                    QP(n,m)=1;
                    if nargout > 1
                        if i==j && j==k && k==l
                            QPinv(m,n)=1;
                        elseif i==j && j==k && k>l
                            QPinv(m,n)=1/4;
                        elseif i>j && j==k && k==l
                            QPinv(m,n)=1/4;
                        elseif i==j && j>k && k==l
                            QPinv(m,n) = 1/6;
                        elseif i>j && j>k && k==l
                            QPinv(m,n) = 1/12;
                        elseif i>j && j==k && k>l
                            QPinv(m,n) = 1/12;
                        elseif i==j && j>k && k>l
                            QPinv(m,n) = 1/12;
                        elseif i>j && j>k && k>l
                            QPinv(m,n) = 1/24;                    
                        end
                    end
                end
                counti = counti+1;
            end
        end
    end
end
%QPinv = (transpose(QP)*QP)\transpose(QP);

function m = mue(p,i,j,k,l)
    % Auxiliary expression, see page 118 of Meijer (2005)
     m = i + (j-1)*p + 1/2*(k-1)*p^2 + 1/6*(l-1)*p^3 - 1/2*j*(j-1) + 1/6*k*(k-1)*(k-2) - 1/24*l*(l-1)*(l-2)*(l-3) - 1/2*(k-1)^2*p + 1/6*(l-1)^3*p - 1/4*(l-1)*(l-2)*p^2 - 1/4*l*(l-1)*p + 1/6*(l-1)*p;
     m = round(m);
end


end