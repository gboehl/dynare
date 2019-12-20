%
% bivmom.m		Date: 1/11/2004
% This Matlab program computes the product moment of X_1^{p_1}X_2^{p_2},
% where X_1 and X_2 are standard bivariate normally distributed.
% n : dimension of X
% rho: correlation coefficient between X_1 and X_2
% Reference: Kotz, Balakrishnan, and Johnson (2000), Continuous Multivariate
%            Distributions, Vol. 1, p.261
% Note that there is a typo in Eq.(46.25), there should be an extra rho in front 
% of the equation.
% Usage: bivmom(p,rho)
%
function [y,dy] = bivmom(p,rho)
s1 = p(1);
s2 = p(2);
rho2 = rho^2;
if nargout > 1
    drho2 = 2*rho;
end
if rem(s1+s2,2)==1
   y = 0;
   return
end
r = fix(s1/2);
s = fix(s2/2);
y = 1;
c = 1;
if nargout > 1
    dy = 0;
    dc = 0;
end
odd = 2*rem(s1,2);
for j=1:min(r,s)
    if nargout > 1
        dc = 2*dc*(r+1-j)*(s+1-j)*rho2/(j*(2*j-1+odd)) + 2*c*(r+1-j)*(s+1-j)*drho2/(j*(2*j-1+odd));    
    end
    c = 2*c*(r+1-j)*(s+1-j)*rho2/(j*(2*j-1+odd));
    y = y+c;
    if nargout > 1
        dy = dy + dc;
    end
end
if odd
   if nargout > 1
    dy = y + dy*rho;
   end
   y = y*rho;
end
y = prod([1:2:s1])*prod([1:2:s2])*y;
if nargout > 1
    dy = prod([1:2:s1])*prod([1:2:s2])*dy;
end

