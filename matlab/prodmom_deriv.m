function [y,dy] = prodmom_deriv(V,ii,nu,dV,dC)
% Computes the product moments (and its derivatives with respect to standard
% errors and correlation parameters) of normally distributed variables, i.e.
% this function computes the product moment of 
% X_{i_1}^{nu_1}X_{i_2}^{nu_2}...X_{i_m}^{nu_m}, where X_{i_j} are elements
% from X ~ N(0_n,V) and V is positive semidefinite.
% Example: To get E[X_2X_4^3X_7^2], use prodmom_deriv(V,[2 4 7],[1 3 2])
% =========================================================================
% INPUTS
%  V    [n by n]     covariance matrix of X (needs to be positive semidefinite)
%  ii   [m by 1]     vector of i_j
%  nu   [nu_m by 1]  power of X_{i_j}
%  dV   [n by n by stderrparam_nbr+corrparam_nbr] derivative of V with respect
%                                                 to selected standard error (stderr)
%                                                 and correlation (corr) parameters
%  dC   [n by n by stderrparam_nbr+corrparam_nbr] derivative of Correlation matrix C with respect
%                                                 to selected standard error (stderr)
%                                                 and correlation (corr) parameters
% -------------------------------------------------------------------------
% OUTPUTS
%  y    [1 by 1]    product moment E[X_{i_1}^{nu_1}X_{i_2}^{nu_2}...X_{i_m}^{nu_m}]
%  dy   [1 by stderrparam_nbr+corrparam_nbr] derivatives of y wrt to selected
%                                            standard error and corr parameters
% -------------------------------------------------------------------------
% This function is based upon prodmom.m which is part of replication codes
% of the following paper:
% Kan, R.: "From moments of sum to moments of product." Journal of 
% Multivariate Analysis, 2008, vol. 99, issue 3, pages 542-554.
% prodmom.m can be retrieved from http://www-2.rotman.utoronto.ca/~kan/papers/prodmom.zip
% Further references:
%  Triantafyllopoulos (2003) On the Central Moments of the Multidimensional
%  Gaussian Distribution, Mathematical Scientist
%  Kotz, Balakrishnan, and Johnson (2000), Continuous Multivariate 
%  Distributions, Vol. 1, p.261
% =========================================================================
% Copyright (C) 2008-2015 Raymond Kan <kan@chass.utoronto.ca>
% Copyright (C) 2019-2020 Dynare Team
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
if nargin<3
    nu = ones(size(ii));
end
s = sum(nu);
if s==0
    y = 1;
    if nargout > 1
        dy = zeros(1,1,size(dV,3));
    end
    return
end
if rem(s,2)==1
    y  = 0;
    if nargout > 1
        dy = zeros(1,1,size(dV,3));
    end
    return
end
nuz = nu==0;
nu(nuz) = [];
ii(nuz) = [];
m = length(ii);
V = V(ii,ii);
if nargout > 1
    dV = dV(ii,ii,:);
end
s2 = s/2;
%
%  Use univariate normal results
%
if m==1
    y = V^s2*prod([1:2:s-1]);
    if nargout > 1
        dy = s2*V^(s2-1)*dV*prod([1:2:s-1]);
        dy = reshape(dy,1,size(dV,3));
    end
   return
end
%
%  Use bivariate normal results when there are only two distinct indices
%
if m==2
    if V(1,1)==0 || V(2,2)==0
        y=0;
        if nargout>1
            dy=zeros(1,size(dV,3));
        end
        return
    end
    rho = V(1,2)/sqrt(V(1,1)*V(2,2));
    if nargout > 1
        drho = dC(ii(1),ii(2),:);
        [tmp,dtmp] = bivmom(nu,rho);
        dy = (nu(1)/2)*V(1,1)^(nu(1)/2-1)*dV(1,1,:) * V(2,2)^(nu(2)/2) * tmp...
           + V(1,1)^(nu(1)/2) * (nu(2)/2)*V(2,2)^(nu(2)/2-1)*dV(2,2,:) * tmp...
           + V(1,1)^(nu(1)/2) * V(2,2)^(nu(2)/2) * dtmp * drho;
        dy = reshape(dy,1,size(dV,3));
    else
        tmp = bivmom(nu,rho);
    end
    y = V(1,1)^(nu(1)/2)*V(2,2)^(nu(2)/2)*tmp;
    return  
end
%
%  Regular case
%
[nu,inu] = sort(nu,2,'descend');
V = V(inu,inu);          % Extract only the relevant part of V
x = zeros(1,m);
V = V./2;
nu2 = nu./2;
p = 2;
q = nu2*V*nu2';
y = 0;
if nargout > 1
    dV = dV(inu,inu,:);      % Extract only the relevant part of dV
    dV = dV./2;
    %dq = nu2*dV*nu2';
    %dq = multiprod(multiprod(nu2,dV),nu2');
    dq = NaN(size(q,1), size(q,2), size(dV,3));
    for jp = 1:size(dV,3)
        dq(:,:,jp) = nu2*dV(:,:,jp)*nu2';
    end
    dy = 0;
end
for i=1:fix(prod(nu+1)/2)
    y = y+p*q^s2;
    if nargout > 1
        dy = dy+p*s2*q^(s2-1)*dq;
    end
    for j=1:m
        if x(j)<nu(j)
            x(j) = x(j)+1;
            p = -round(p*(nu(j)+1-x(j))/x(j));
            q = q-2*(nu2-x)*V(:,j)-V(j,j);
            if nargout > 1
                %dq = dq-2*(nu2-x)*dV(:,j,:)-dV(j,j,:);
                %dq = dq-2*multiprod((nu2-x),dV(:,j,:))-dV(j,j,:);
                for jp=1:size(dV,3)
                    dq(:,:,jp) = dq(:,:,jp)-2*(nu2-x)*dV(:,j,jp)-dV(j,j,jp);
                end
            end
            break
        else
            x(j) = 0;
            if rem(nu(j),2)==1
                p = -p;
            end
            if nargout > 1
                %dq = dq+2*nu(j)*multiprod((nu2-x),dV(:,j,:))-nu(j)^2*dV(j,j,:);
                for jp=1:size(dV,3)
                    dq(:,:,jp) = dq(:,:,jp)+2*nu(j)*(nu2-x)*dV(:,j,jp)-nu(j)^2*dV(j,j,jp);
                end
            end
            q = q+2*nu(j)*(nu2-x)*V(:,j)-nu(j)^2*V(j,j);
        end
    end
end
y = y/prod([1:s2]);
if nargout > 1
    dy = dy/prod([1:s2]);
    dy = reshape(dy,1,size(dV,3));
end