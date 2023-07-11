function p = draw(o)

% Return a random draw from the prior distribution.
%
% INPUTS
% - o    [dprior]
%
% OUTPUTS
% - p    [double]   m×1 vector, random draw from the prior distribution (m is the number of estimated parameters).
%
% REMARKS
% None.
%
% EXAMPLE
%
% >> Prior = dprior(bayestopt_, options_.prior_trunc);
% >> d = Prior.draw()

% Copyright © 2023 Dynare Team
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

p = NaN(rows(o.lb), 1);
if o.isuniform
    p(o.iduniform) = rand(length(o.iduniform),1).*(o.p4(o.iduniform)-o.p3(o.iduniform)) + o.p3(o.iduniform);
    oob = find( (p(o.iduniform)>o.ub(o.iduniform)) | (p(o.iduniform)<o.lb(o.iduniform)));
    while ~isempty(oob)
        p(o.iduniform) = rand(length(o.iduniform), 1).*(o.p4(o.iduniform)-o.p3(o.iduniform)) + o.p3(o.iduniform);
        oob = find( (p(o.iduniform)>o.ub(o.iduniform)) | (p(o.iduniform)<o.lb(o.iduniform)));
    end
end
if o.isgaussian
    p(o.idgaussian) = randn(length(o.idgaussian), 1).*o.p7(o.idgaussian) + o.p6(o.idgaussian);
    oob = find( (p(o.idgaussian)>o.ub(o.idgaussian)) | (p(o.idgaussian)<o.lb(o.idgaussian)));
    while ~isempty(oob)
        p(o.idgaussian(oob)) = randn(length(o.idgaussian(oob)), 1).*o.p7(o.idgaussian(oob)) + o.p6(o.idgaussian(oob));
        oob = find( (p(o.idgaussian)>o.ub(o.idgaussian)) | (p(o.idgaussian)<o.lb(o.idgaussian)));
    end
end
if o.isgamma
    p(o.idgamma) = gamrnd(o.p6(o.idgamma), o.p7(o.idgamma))+o.p3(o.idgamma);
    oob = find( (p(o.idgamma)>o.ub(o.idgamma)) | (p(o.idgamma)<o.lb(o.idgamma)));
    while ~isempty(oob)
        p(o.idgamma(oob)) = gamrnd(o.p6(o.idgamma(oob)), o.p7(o.idgamma(oob)))+o.p3(o.idgamma(oob));
        oob = find( (p(o.idgamma)>o.ub(o.idgamma)) | (p(o.idgamma)<o.lb(o.idgamma)));
    end
end
if o.isbeta
    p(o.idbeta) = (o.p4(o.idbeta)-o.p3(o.idbeta)).*betarnd(o.p6(o.idbeta), o.p7(o.idbeta))+o.p3(o.idbeta);
    oob = find( (p(o.idbeta)>o.ub(o.idbeta)) | (p(o.idbeta)<o.lb(o.idbeta)));
    while ~isempty(oob)
        p(o.idbeta(oob)) = (o.p4(o.idbeta(oob))-o.p3(o.idbeta(oob))).*betarnd(o.p6(o.idbeta(oob)), o.p7(o.idbeta(oob)))+o.p3(o.idbeta(oob));
        oob = find( (p(o.idbeta)>o.ub(o.idbeta)) | (p(o.idbeta)<o.lb(o.idbeta)));
    end
end
if o.isinvgamma1
    p(o.idinvgamma1) = ...
        sqrt(1./gamrnd(o.p7(o.idinvgamma1)/2, 2./o.p6(o.idinvgamma1)))+o.p3(o.idinvgamma1);
    oob = find( (p(o.idinvgamma1)>o.ub(o.idinvgamma1)) | (p(o.idinvgamma1)<o.lb(o.idinvgamma1)));
    while ~isempty(oob)
        p(o.idinvgamma1(oob)) = ...
            sqrt(1./gamrnd(o.p7(o.idinvgamma1(oob))/2, 2./o.p6(o.idinvgamma1(oob))))+o.p3(o.idinvgamma1(oob));
        oob = find( (p(o.idinvgamma1)>o.ub(o.idinvgamma1)) | (p(o.idinvgamma1)<o.lb(o.idinvgamma1)));
    end
end
if o.isinvgamma2
    p(o.idinvgamma2) = ...
        1./gamrnd(o.p7(o.idinvgamma2)/2, 2./o.p6(o.idinvgamma2))+o.p3(o.idinvgamma2);
    oob = find( (p(o.idinvgamma2)>o.ub(o.idinvgamma2)) | (p(o.idinvgamma2)<o.lb(o.idinvgamma2)));
    while ~isempty(oob)
        p(o.idinvgamma2(oob)) = ...
            1./gamrnd(o.p7(o.idinvgamma2(oob))/2, 2./o.p6(o.idinvgamma2(oob)))+o.p3(o.idinvgamma2(oob));
        oob = find( (p(o.idinvgamma2)>o.ub(o.idinvgamma2)) | (p(o.idinvgamma2)<o.lb(o.idinvgamma2)));
    end
end
if o.isweibull
    p(o.idweibull) = wblrnd(o.p7(o.idweibull), o.p6(o.idweibull)) + o.p3(o.idweibull);
    oob = find( (p(o.idweibull)>o.ub(o.idweibull)) | (p(o.idweibull)<o.lb(o.idweibull)));
    while ~isempty(oob)
        p(o.idweibull(oob)) = wblrnd(o.p7(o.idweibull(oob)), o.p6(o.idweibull(oob)))+o.p3(o.idweibull(oob));
        oob = find( (p(o.idweibull)>o.ub(o.idweibull)) | (p(o.idweibull)<o.lb(o.idweibull)));
    end
end

return % --*-- Unit tests --*--

%@test:1
% Fill global structures with required fields...
prior_trunc = 1e-10;
p0 = repmat([1; 2; 3; 4; 5; 6; 8], 2, 1);    % Prior shape
p1 = .4*ones(14,1);                          % Prior mean
p2 = .2*ones(14,1);                          % Prior std.
p3 = NaN(14,1);
p4 = NaN(14,1);
p5 = NaN(14,1);
p6 = NaN(14,1);
p7 = NaN(14,1);

for i=1:14
   switch p0(i)
     case 1
       % Beta distribution
       p3(i) = 0;
       p4(i) = 1;
       [p6(i), p7(i)] = beta_specification(p1(i), p2(i)^2, p3(i), p4(i));
       p5(i) = compute_prior_mode([p6(i) p7(i)], 1);
     case 2
       % Gamma distribution
       p3(i) = 0;
       p4(i) = Inf;
       [p6(i), p7(i)] = gamma_specification(p1(i), p2(i)^2, p3(i), p4(i));
       p5(i) = compute_prior_mode([p6(i) p7(i)], 2);
     case 3
       % Normal distribution
       p3(i) = -Inf;
       p4(i) = Inf;
       p6(i) = p1(i);
       p7(i) = p2(i);
       p5(i) = p1(i);
     case 4
       % Inverse Gamma (type I) distribution
       p3(i) = 0;
       p4(i) = Inf;
       [p6(i), p7(i)] = inverse_gamma_specification(p1(i), p2(i)^2, p3(i), 1, false);
       p5(i) = compute_prior_mode([p6(i) p7(i)], 4);
     case 5
       % Uniform distribution
       [p1(i), p2(i), p6(i), p7(i)] = uniform_specification(p1(i), p2(i), p3(i), p4(i));
       p3(i) = p6(i);
       p4(i) = p7(i);
       p5(i) = compute_prior_mode([p6(i) p7(i)], 5);
     case 6
       % Inverse Gamma (type II) distribution
       p3(i) = 0;
       p4(i) = Inf;
       [p6(i), p7(i)] = inverse_gamma_specification(p1(i), p2(i)^2, p3(i), 2, false);
       p5(i) = compute_prior_mode([p6(i) p7(i)], 6);
     case 8
       % Weibull distribution
       p3(i) = 0;
       p4(i) = Inf;
       [p6(i), p7(i)] = weibull_specification(p1(i), p2(i)^2, p3(i));
       p5(i) = compute_prior_mode([p6(i) p7(i)], 8);
     otherwise
       error('This density is not implemented!')
   end
end

BayesInfo.pshape = p0;
BayesInfo.p1 = p1;
BayesInfo.p2 = p2;
BayesInfo.p3 = p3;
BayesInfo.p4 = p4;
BayesInfo.p5 = p5;
BayesInfo.p6 = p6;
BayesInfo.p7 = p7;

ndraws = 1e5;
m0 = BayesInfo.p1; %zeros(14,1);
v0 = diag(BayesInfo.p2.^2); %zeros(14);

% Call the tested routine
try
   % Instantiate dprior object
   o = dprior(BayesInfo, prior_trunc, false);
   % Do simulations in a loop and estimate recursively the mean and the variance.
   for i = 1:ndraws
        draw = o.draw();
        m1 = m0 + (draw-m0)/i;
        m2 = m1*m1';
        v0 = v0 + ((draw*draw'-m2-v0) + (i-1)*(m0*m0'-m2'))/i;
        m0 = m1;
   end
   t(1) = true;
catch
    t(1) = false;
end

if t(1)
    t(2) = all(abs(m0-BayesInfo.p1)<3e-3);
    t(3) = all(all(abs(v0-diag(BayesInfo.p2.^2))<5e-3));
end
T = all(t);
%@eof:1
