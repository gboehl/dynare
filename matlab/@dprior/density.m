function [lpd, dlpd, d2lpd, info] = density(o, x)

% Evaluate the logged prior density at x.
%
% INPUTS
% - o       [dprior]
% - x       [double]   m×1 vector, point where the prior density is evaluated.
%
% OUTPUTS
% - lpd     [double]   scalar, value of the logged prior density at x.
% - dlpd    [double]   m×1 vector, first order derivatives.
% - d2lpd   [double]   m×1 vector, second order derivatives.
%
% REMARKS
% Second order derivatives holder, d2lpd, has the same rank and shape than dlpd because the priors are
% independent (we would have to use a matrix if non orthogonal priors were allowed in Dynare).
%
% EXAMPLE
%
% >> Prior = dprior(bayestopt_, options_.prior_trunc);
% >> lpd = Prior.dsensity(x)

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

lpd = 0.0;
if nargout>1
    dlpd = zeros(1, length(x));
    if nargout>2
        d2lpd = dlpd;
        if nargout>3
            info = [];
        end
    end
end
if o.isuniform
    if any(x(o.iduniform)-o.p3(o.iduniform)<0) || any(x(o.iduniform)-o.p4(o.iduniform)>0)
        lpd = -Inf ;
        if nargout==4
            info = o.iduniform((x(o.iduniform)-o.p3(o.iduniform)<0) || (x(o.iduniform)-o.p4(o.iduniform)>0));
        end
        return
    end
    lpd = lpd - sum(log(o.p4(o.iduniform)-o.p3(o.iduniform))) ;
    if nargout>1
        dlpd(o.iduniform) = zeros(length(o.iduniform), 1);
        if nargout>2
            d2lpd(o.iduniform) = zeros(length(o.iduniform), 1);
        end
    end
end
if o.isgaussian
    switch nargout
      case 1
        lpd = lpd + sum(lpdfnorm(x(o.idgaussian), o.p6(o.idgaussian), o.p7(o.idgaussian)));
      case 2
        [tmp, dlpd(o.idgaussian)] = lpdfnorm(x(o.idgaussian), o.p6(o.idgaussian), o.p7(o.idgaussian));
        lpd = lpd + sum(tmp);
      case {3,4}
        [tmp, dlpd(o.idgaussian), d2lpd(o.idgaussian)] = lpdfnorm(x(o.idgaussian), o.p6(o.idgaussian), o.p7(o.idgaussian));
        lpd = lpd + sum(tmp);
    end
end
if o.isgamma
    switch nargout
      case 1
        lpd = lpd + sum(lpdfgam(x(o.idgamma)-o.p3(o.idgamma), o.p6(o.idgamma), o.p7(o.idgamma)));
        if isinf(lpd), return, end
      case 2
        [tmp, dlpd(o.idgamma)] = lpdfgam(x(o.idgamma)-o.p3(o.idgamma), o.p6(o.idgamma), o.p7(o.idgamma));
        lpd = lpd + sum(tmp);
        if isinf(lpd), return, end
      case 3
        [tmp, dlpd(o.idgamma), d2lpd(o.idgamma)] = lpdfgam(x(o.idgamma)-o.p3(o.idgamma), o.p6(o.idgamma), o.p7(o.idgamma));
        lpd = lpd + sum(tmp);
        if isinf(lpd), return, end
      case 4
        [tmp, dlpd(o.idgamma), d2lpd(o.idgamma)] = lpdfgam(x(o.idgamma)-o.p3(o.idgamma), o.p6(o.idgamma), o.p7(o.idgamma));
        lpd = lpd + sum(tmp);
        if isinf(lpd)
            info = o.idgamma(isinf(tmp));
            return
        end
    end
end
if o.isbeta
    switch nargout
      case 1
        lpd = lpd + sum(lpdfgbeta(x(o.idbeta), o.p6(o.idbeta), o.p7(o.idbeta), o.p3(o.idbeta), o.p4(o.idbeta)));
        if isinf(lpd), return, end
      case 2
        [tmp, dlpd(o.idbeta)] = lpdfgbeta(x(o.idbeta), o.p6(o.idbeta), o.p7(o.idbeta), o.p3(o.idbeta), o.p4(o.idbeta));
        lpd = lpd + sum(tmp);
        if isinf(lpd), return, end
      case 3
        [tmp, dlpd(o.idbeta), d2lpd(o.idbeta)] = lpdfgbeta(x(o.idbeta), o.p6(o.idbeta), o.p7(o.idbeta), o.p3(o.idbeta), o.p4(o.idbeta));
        lpd = lpd + sum(tmp);
        if isinf(lpd), return, end
      case 4
        [tmp, dlpd(o.idbeta), d2lpd(o.idbeta)] = lpdfgbeta(x(o.idbeta), o.p6(o.idbeta), o.p7(o.idbeta), o.p3(o.idbeta), o.p4(o.idbeta));
        lpd = lpd + sum(tmp);
        if isinf(lpd)
            info = o.idbeta(isinf(tmp));
            return
        end
    end
end
if o.isinvgamma1
    switch nargout
      case 1
        lpd = lpd + sum(lpdfig1(x(o.idinvgamma1)-o.p3(o.idinvgamma1), o.p6(o.idinvgamma1), o.p7(o.idinvgamma1)));
        if isinf(lpd), return, end
      case 2
        [tmp, dlpd(o.idinvgamma1)] = lpdfig1(x(o.idinvgamma1)-o.p3(o.idinvgamma1), o.p6(o.idinvgamma1), o.p7(o.idinvgamma1));
        lpd = lpd + sum(tmp);
        if isinf(lpd), return, end
      case 3
        [tmp, dlpd(o.idinvgamma1), d2lpd(o.idinvgamma1)] = lpdfig1(x(o.idinvgamma1)-o.p3(o.idinvgamma1), o.p6(o.idinvgamma1), o.p7(o.idinvgamma1));
        lpd = lpd + sum(tmp);
        if isinf(lpd), return, end
      case 4
        [tmp, dlpd(o.idinvgamma1), d2lpd(o.idinvgamma1)] = lpdfig1(x(o.idinvgamma1)-o.p3(o.idinvgamma1), o.p6(o.idinvgamma1), o.p7(o.idinvgamma1));
        lpd = lpd + sum(tmp);
        if isinf(lpd)
            info = o.idinvgamma1(isinf(tmp));
            return
        end
    end
end
if o.isinvgamma2
    switch nargout
      case 1
        lpd = lpd + sum(lpdfig2(x(o.idinvgamma2)-o.p3(o.idinvgamma2), o.p6(o.idinvgamma2), o.p7(o.idinvgamma2)));
        if isinf(lpd), return, end
      case 2
        [tmp, dlpd(o.idinvgamma2)] = lpdfig2(x(o.idinvgamma2)-o.p3(o.idinvgamma2), o.p6(o.idinvgamma2), o.p7(o.idinvgamma2));
        lpd = lpd + sum(tmp);
        if isinf(lpd), return, end
      case 3
        [tmp, dlpd(o.idinvgamma2), d2lpd(o.idinvgamma2)] = lpdfig2(x(o.idinvgamma2)-o.p3(o.idinvgamma2), o.p6(o.idinvgamma2), o.p7(o.idinvgamma2));
        lpd = lpd + sum(tmp);
        if isinf(lpd), return, end
      case 4
        [tmp, dlpd(o.idinvgamma2), d2lpd(o.idinvgamma2)] = lpdfig2(x(o.idinvgamma2)-o.p3(o.idinvgamma2), o.p6(o.idinvgamma2), o.p7(o.idinvgamma2));
        lpd = lpd + sum(tmp);
        if isinf(lpd)
            info = o.idinvgamma2(isinf(tmp));
            return
        end
    end
end
if o.isweibull
    switch nargout
      case 1
        lpd = lpd + sum(lpdfgweibull(x(o.idweibull), o.p6(o.idweibull), o.p7(o.idweibull)));
        if isinf(lpd), return, end
      case 2
        [tmp, dlpd(o.idweibull)] = lpdfgweibull(x(o.idweibull), o.p6(o.idweibull), o.p7(o.idweibull));
        lpd = lpd + sum(tmp);
        if isinf(lpd), return, end
      case 3
        [tmp, dlpd(o.idweibull), d2lpd(o.idweibull)] = lpdfgweibull(x(o.idweibull), o.p6(o.idweibull), o.p7(o.idweibull));
        lpd = lpd + sum(tmp);
        if isinf(lpd), return, end
      case 4
        [tmp, dlpd(o.idweibull), d2lpd(o.idweibull)] = lpdfgweibull(x(o.idweibull), o.p6(o.idweibull), o.p7(o.idweibull));
        lpd = lpd + sum(tmp);
        if isinf(lpd)
            info = o.idweibull(isinf(tmp));
            return
        end
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

% Call the tested routine
try
    Prior = dprior(BayesInfo, prior_trunc, false);

    % Compute density at the prior mode
    lpdstar = Prior.density(p5);

    % Draw random deviates in a loop and evaluate the density.
    LPD = NaN(10000,1);
    parfor i = 1:10000
        x = Prior.draw();
        LPD(i) = Prior.density(x);
    end
    t(1) = true;
catch
    t(1) = false;
end

if t(1)
    t(2) = all(LPD<=lpdstar);
end
T = all(t);
%@eof:1

%@test:2
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

% Call the tested routine
try
    Prior = dprior(BayesInfo, prior_trunc, false);
    mu = NaN(14,1);
    std = NaN(14,1);

    for i=1:14
        % Define conditional density (it's also a marginal since priors are orthogonal)
        f = @(x) exp(Prior.densities(substitute(p5, i, x)));
        % TODO: Check the version of Octave we use (integral is available as a compatibility wrapper in latest Octave version)
        m = integral(f, p3(i), p4(i));
        density = @(x) f(x)/m; % rescaling is required since the probability mass depends on the conditioning.
        % Compute the conditional expectation
        mu(i) = integral(@(x) x.*density(x), p3(i), p4(i));
        std(i) = sqrt(integral(@(x) ((x-mu(i)).^2).*density(x), p3(i), p4(i)));
    end

    t(1) = true;
catch
    t(1) = false;
end

if t(1)
    t(2) = all(abs(mu-.4)<1e-6);
    t(3) = all(abs(std-.2)<1e-6);
end
T = all(t);
%@eof:2
