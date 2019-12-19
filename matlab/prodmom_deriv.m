%
% prodmom.m		Date: 4/29/2006
% This Matlab program computes the product moment of X_{i_1}^{nu_1}X_{i_2}^{nu_2}...X_{i_m}^{nu_m},
% where X_{i_j} are elements from X ~ N(0_n,V).  
% V only needs to be positive semidefinite.
% V: variance-covariance matrix of X
% ii: vector of i_j
% nu: power of X_{i_j} 
% Reference: Triantafyllopoulos (2003) On the Central Moments of the Multidimensional
%            Gaussian Distribution, Mathematical Scientist
%            Kotz, Balakrishnan, and Johnson (2000), Continuous Multivariate
%            Distributions, Vol. 1, p.261
% Note that there is a typo in Eq.(46.25), there should be an extra rho in front 
% of the equation.
% Usage: prodmom(V,[i1 i2 ... ir],[nu1 nu2 ... nur])
% Example: To get E[X_2X_4^3X_7^2], use prodmom(V,[2 4 7],[1 3 2])
%
function dy = prodmom_deriv(V,ii,nu,dV,dC);
if nargin<3
   nu = ones(size(ii));
end
s = sum(nu);
if s==0
   dy = zeros(1,1,size(dV,3));
   return
end
if rem(s,2)==1
    dy = zeros(1,1,size(dV,3));
   return
end
nuz = nu==0;
nu(nuz) = [];
ii(nuz) = [];
m = length(ii);
V = V(ii,ii);
dV = dV(ii,ii,:);
s2 = s/2;
%
%  Use univariate normal results
%
if m==1
   dy = s2*V^(s2-1)*dV*prod([1:2:s-1]);
   dy = reshape(dy,1,size(dV,3));
   return
end
%
%  Use bivariate normal results when there are only two distinct indices
%
if m==2
   rho = V(1,2)/sqrt(V(1,1)*V(2,2));
   drho = dC(ii(1),ii(2),:);
   [tmp,dtmp] = bivmom(nu,rho);
   dy =   (nu(1)/2)*V(1,1)^(nu(1)/2-1)*dV(1,1,:) * V(2,2)^(nu(2)/2) * tmp...
       + V(1,1)^(nu(1)/2) * (nu(2)/2)*V(2,2)^(nu(2)/2-1)*dV(2,2,:) * tmp...
       + V(1,1)^(nu(1)/2) * V(2,2)^(nu(2)/2) * dtmp * drho;
   dy = reshape(dy,1,size(dV,3));
   return  
end
%
%  Regular case
%
[nu,inu] = sort(nu,2,'descend');
V = V(inu,inu);          % Extract only the relevant part of V
dV = dV(inu,inu,:);          % Extract only the relevant part of dV
x = zeros(1,m);
V = V./2;
dV = dV./2;
nu2 = nu./2;
p = 2;
q = nu2*V*nu2';
%dq = nu2*dV*nu2';
%dq = multiprod(multiprod(nu2,dV),nu2');
dq = NaN(size(q,1), size(q,2), size(dV,3));
for jp = 1:size(dV,3)
    dq(:,:,jp) = nu2*dV(:,:,jp)*nu2';
end
dy = 0;
for i=1:fix(prod(nu+1)/2)
    dy = dy+p*s2*q^(s2-1)*dq;
    for j=1:m
        if x(j)<nu(j)
           x(j) = x(j)+1;
           p = -round(p*(nu(j)+1-x(j))/x(j));           
           %dq = dq-2*(nu2-x)*dV(:,j,:)-dV(j,j,:);
           %dq = dq-2*multiprod((nu2-x),dV(:,j,:))-dV(j,j,:);
           for jp=1:size(dV,3)
               dq(:,:,jp) = dq(:,:,jp)-2*(nu2-x)*dV(:,j,jp)-dV(j,j,jp);
           end
           q = q-2*(nu2-x)*V(:,j)-V(j,j);
           break
        else
           x(j) = 0;
           if rem(nu(j),2)==1
              p = -p;
           end           
           %dq = dq+2*nu(j)*multiprod((nu2-x),dV(:,j,:))-nu(j)^2*dV(j,j,:);
           for jp=1:size(dV,3)
               dq(:,:,jp) = dq(:,:,jp)+2*nu(j)*(nu2-x)*dV(:,j,jp)-nu(j)^2*dV(j,j,jp);
           end
           q = q+2*nu(j)*(nu2-x)*V(:,j)-nu(j)^2*V(j,j);
        end
    end
end
dy = dy/prod([1:s2]);
dy = reshape(dy,1,size(dV,3));