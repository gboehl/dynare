function [DP6,DP6inv] = Q6_plication(p)
% Computes the 6-way duplication Matrix DP6 (and its Moore-Penrose inverse)
% such that for any p-dimensional vector x:
% y=kron(kron(kron(kron(kron(x,x),x,x),x),x)=DP6*z
% where z is of dimension np=p*(p+1)*(p+2)*(p+3)*(p+4)*(p+5)/(1*2*3*4*5*6)
% and is obtained from y by removing each second and later occurence of the
% same element. This is a generalization of the Duplication matrix.
% Reference: Meijer (2005) - Matrix algebra for higher order moments.
%            Linear Algebra and its Applications, 410,pp. 112-134
% =========================================================================
% INPUTS
%    * p        [integer]    size of vector
% -------------------------------------------------------------------------
% OUTPUTS
%    * DP6      [p^6 by np]  6-way duplication matrix
%    * DP6inv   [np by np]   Moore-Penrose inverse of DP6
% -------------------------------------------------------------------------
% This function is called by
%   * pruned_state_space_system.m
% -------------------------------------------------------------------------
% This function calls
%   * binom_coef (embedded)
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
np = p*(p+1)*(p+2)*(p+3)*(p+4)*(p+5)/(1*2*3*4*5*6);
DP6 = spalloc(p^6,np,p^6);
counti=1;
for i1=1:p
    for i2=i1:p
        for i3=i2:p
            for i4=i3:p
                for i5=i4:p
                    for i6=i5:p
                        idx = pruned_SS.uperm([i6 i5 i4 i3 i2 i1]);
                        for r = 1:size(idx,1)
                            ii1 = idx(r,1); ii2= idx(r,2); ii3=idx(r,3); ii4=idx(r,4); ii5=idx(r,5); ii6=idx(r,6);
                            n = ii1 + (ii2-1)*p + (ii3-1)*p^2 + (ii4-1)*p^3  + (ii5-1)*p^4 + (ii6-1)*p^5;
                            m = mue(p,i6,i5,i4,i3,i2,i1);
                            DP6(n,m)=1;
                        end
                        counti = counti+1;
                    end                
                end
            end
        end
    end
end
if nargout==2
DP6inv = (transpose(DP6)*DP6)\transpose(DP6);
end

function m = mue(p,i1,i2,i3,i4,i5,i6)
% Auxiliary expression, see page 122 of Meijer (2005)
     m = binom_coef(p,6,1) - binom_coef(p,1,i1+1) - binom_coef(p,2,i2+1) - binom_coef(p,3,i3+1) - binom_coef(p,4,i4+1) - binom_coef(p,5,i5+1) - binom_coef(p,6,i6+1);
     m = round(m);
end

function N = binom_coef(p,q,i)
% Auxiliary expression for binomial coefficients, see page 119 of Meijer (2005)
    t = q; r =p+q-i;
    if t==0
        N=1;
    else
        N=1;
        for h = 0:(t-1)
            N = N*(r-h);
        end
        N=N/factorial(t);
    end
end
end