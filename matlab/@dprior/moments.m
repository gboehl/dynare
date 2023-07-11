function o = moments(o, name)

% Compute the prior moments.
%
% INPUTS
% - o       [dprior]
%
% OUTPUTS
% - o       [dprior]

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

switch name
  case 'mean'
    m = o.p1;
  case 'median'
    m = o.p11;
  case 'std'
    m = o.p2;
  case 'mode'
    m = o.p5;
  otherwise
    error('%s is not an implemented moemnt.', name)
end
id = isnan(m);
if any(id)
    % For some parameters the prior mean is not defined.
    % We compute the first order moment from the
    % hyperparameters, if the hyperparameters are not
    % available an error is thrown.
    if o.isuniform
        jd = intersect(o.iduniform, find(id));
        if ~isempty(jd)
            if any(isnan(o.p3(jd))) || any(isnan(o.p4(jd)))
                error('dprior::mean: Some hyperparameters are missing (uniform distribution).')
            end
            switch name
              case 'mean'
                m(jd) = o.p3(jd) + .5*(o.p4(jd)-o.p3(jd));
              case 'median'
                m(jd) = o.p3(jd) + .5*(o.p4(jd)-o.p3(jd));
              case 'std'
                m(jd) = (o.p4(jd)-o.p3(jd))/sqrt(12);
              case 'mode' % Actually we have a continuum of modes with the uniform distribution.
                m(jd) = o.p3(jd) + .5*(o.p4(jd)-o.p3(jd));
            end
        end
    end
    if o.isgaussian
        jd = intersect(o.idgaussian, find(id));
        if ~isempty(jd)
            if any(isnan(o.p6(jd))) || any(isnan(o.p7(jd)))
                error('dprior::mean: Some hyperparameters are missing (gaussian distribution).')
            end
            switch name
              case 'mean'
                m(jd) = o.p6(jd);
              case 'median'
                m(jd) = o.p6(jd);
              case 'std'
                m(jd) = o.p7(jd);
              case 'mode' % Actually we have a continuum of modes with the uniform distribution.
                m(jd) = o.p6(jd);
            end
        end
    end
    if o.isgamma
        jd = intersect(o.idgamma, find(id));
        if ~isempty(jd)
            if any(isnan(o.p6(jd))) || any(isnan(o.p7(jd))) || any(isnan(o.p3(jd)))
                error('dprior::mean: Some hyperparameters are missing (gamma distribution).')
            end
            % α → o.p6, β → o.p7
            switch name
              case 'mean'
                m(jd) = o.p3(jd) + o.p6(jd).*o.p7(jd);
              case 'median'
                m(jd) = o.p3(jd) + gaminv(.5, o.p6(jd), o.p7(jda));
              case 'std'
                m(jd) = sqrt(o.p6(jd)).*o.p7(jd);
              case 'mode'
                m(jd) = 0;
                hd = o.p6(jd)>1;
                m(jd(hd)) = (o.p6(jd(hd))-1).*o.p7(jd(hd));
            end
        end
    end
    if o.isbeta
        jd = intersect(o.idbeta, find(id));
        if ~isempty(jd)
            if any(isnan(o.p6(jd))) || any(isnan(o.p7(jd))) || any(isnan(o.p3(jd))) || any(isnan(o.p4(jd)))
                error('dprior::mean: Some hyperparameters are missing (beta distribution).')
            end
            % α → o.p6, β → o.p7
            switch name
              case 'mean'
                m(jd) = o.p3(jd) + (o.p6(jd)./(o.p6(jd)+o.p7(jd))).*(o.p4(jd)-o.p3(jd));
              case 'median'
                m(jd) = o.p3(jd) + betainv(.5, o.p6(jd), o.p7(jd)).*(o.p4(jd)-o.p3(jd));
              case 'std'
                m(jd) = (o.p4(jd)-o.p3(jd)).*sqrt(o.p6(jd).*o.p7(jd)./((o.p6(jd)+o.p7(jd)).^2.*(o.p6(jd)+o.p7(jd)+1)));
              case 'mode'
                h0 = true(jd, 1);
                h1 = o.p6(jd)<=1 & o.p7(jd)>1; h0 = h0 & ~h1;
                h2 = o.p7(jd)<=1 & o.p6(jd)>1; h0 = h0 & ~h2;
                h3 = o.p6(jd)<1 & o.p7(jd)<1; h0 = h0 & ~h3;
                h4 = ismembertol(o.p6(jd), 1) & ismembertol(o.p7(jd),1); h0 = h0 & ~h4;
                m(jd(h1)) = o.p3(jd(h1));                                  % Standard β has a mode at 0
                m(jd(h2)) = o.p4(jd(h2));                                  % Standard β has a mode at 1
                m(jd(h3)) = o.p3(jd(h3));                                  % Standard β is bimodal, we pick the lowest mode (0)
                m(jd(h4)) = o.p3(jd(h4)) + .5*(o.p4(jd(h4))-o.p3(jd(h4))); % Standard β is the uniform distribution (continuum of modes), we pick the mean as the mode
                m(jd(h0)) = o.p3(jd(h0))+(o.p4(jd(h0))-o.p3(jd(h0))).*((o.p6(jd(h0))-1)./(o.p6(jd(h0))+o.p7(jd(h0))-2)); % β distribution is concave and has a unique interior mode.
            end
        end
    end
    if o.isinvgamma1
        jd = intersect(o.idinvgamma1, find(id));
        if ~isempty(jd)
            if any(isnan(o.p6(jd))) || any(isnan(o.p7(jd))) || any(isnan(o.p3(jd)))
                error('dprior::mean: Some hyperparameters are missing (inverse gamma type 1 distribution).')
            end
            % s → o.p6, ν → o.p7
            switch name
              case 'mean'
                m(jd) = o.p3(jd) + sqrt(.5*o.p6(jd)) .*(gamma(.5*(o.p7(jd)-1))./gamma(.5*o.p7(jd)));
              case 'median'
                m(jd) = o.p3(jd) + 1.0/sqrt(gaminv(.5, o.p7(jd)/2.0, 2.0/o.p6(jd)));
              case 'std'
                m(jd) = sqrt( o.p6(jd)./(o.p7(jd)-2)-(.5*o.p6(jd)).*(gamma(.5*(o.p7(jd)-1))./gamma(.5*o.p7(jd))).^2);
              case 'mode'
                m(jd) = sqrt((o.p7(jd)-1)./o.p6(jd));
            end
        end
    end
    if o.isinvgamma2
        jd = intersect(o.idinvgamma2, find(id));
        if ~isempty(jd)
            if any(isnan(o.p6(jd))) || any(isnan(o.p7(jd))) || any(isnan(o.p3(jd)))
                error('dprior::mean: Some hyperparameters are missing (inverse gamma type 2 distribution).')
            end
            % s → o.p6, ν → o.p7
            switch name
              case 'mean'
                m(jd) =  o.p3(jd) + o.p6(jd)./(o.p7(jd)-2);
              case 'median'
                m(jd) = o.p3(jd) + 1.0/gaminv(.5, o.p7(jd)/2.0, 2.0/o.p6(jd));
              case 'std'
                m(jd) = sqrt(2./(o.p7(jd)-4)).*o.p6(jd)./(o.p7(jd)-2);
              case 'mode'
                m(jd) = o.p6(jd)./(o.p7(jd)+2);
            end
        end
    end
    if o.isweibull
        jd = intersect(o.idweibull, find(id));
        if ~isempty(jd)
            if any(isnan(o.p6(jd))) || any(isnan(o.p7(jd))) || any(isnan(o.p3(jd)))
                error('dprior::mean: Some hyperparameters are missing (weibull distribution).')
            end
            % k → o.p6 (shape parameter), λ → o.p7 (scale parameter)
            % See https://en.wikipedia.org/wiki/Weibull_distribution
            switch name
              case 'mean'
                m(jd) =  o.p3(jd) + o.p7(jd).*gamma(1+1./o.p6(jd));
              case 'median'
                m(jd) = o.p3(jd) + o.p7(jd).*log(2).^(1./o.p6(jd));
              case 'std'
                m(jd) = o.p7(jd).*sqrt(gamma(1+2./o.p6(jd))-gamma(1+1./o.p6(jd)).^2);
              case 'mode'
                m(jd) = 0;
                hd = o.p6(jd)>1;
                m(jd(hd)) = o.p3(jd(hd)) + o.p7(jd(hd)).*((o.p6(jd(hd))-1)./o.p6(jd(hd))).^(1./o.p6(jd(hd)));
            end
        end
    end
    switch name
      case 'mean'
        o.p1 = m;
      case 'median'
        o.p11 = m;
      case 'std'
        o.p2 = m;
      case 'mode'
        o.p5 = m;
    end
end

return % --*-- Unit tests --*--

%@test:5
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
    t(1) = true;
catch
    t(1) = false;
end

if t(1)
    t(2) = all(Prior.mean()==.4);
    t(3) = all(ismembertol(Prior.mean(true),.4));
    t(4) = all(ismembertol(Prior.variance(),.04));
    t(5) = all(ismembertol(Prior.variance(true),.04));
end
T = all(t);
%@eof:5
