classdef dprior < handle

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

    properties
        p1 = [];                         % Prior mean.
        p2 = [];                         % Prior stddev.
        p3 = [];                         % Lower bound of the prior support.
        p4 = [];                         % Upper bound of the prior support.
        p5 = [];                         % Prior mode.
        p6 = [];                         % Prior first hyperparameter.
        p7 = [];                         % Prior second hyperparameter.
        p11 = [];                        % Prior median
        lb = [];                         % Truncated prior lower bound.
        ub = [];                         % Truncated prior upper bound.
        name = {};                       % Name of the parameter
        iduniform = [];                  % Index for the uniform priors.
        idgaussian = [];                 % Index for the gaussian priors.
        idgamma = [];                    % Index for the gamma priors.
        idbeta = [];                     % Index for the beta priors.
        idinvgamma1 = [];                % Index for the inverse gamma type 1 priors.
        idinvgamma2 = [];                % Index for the inverse gamma type 2 priors.
        idweibull = [];                  % Index for the weibull priors.
        isuniform = false;
        isgaussian = false;
        isgamma = false;
        isbeta = false;
        isinvgamma1 = false;
        isinvgamma2 = false;
        isweibull = false;
    end

    methods

        function o = dprior(bayestopt_, PriorTrunc, Uniform)
        % Class constructor.
        %
        % INPUTS
        % - bayestopt_    [struct]   Informations about the prior distribution, aka bayestopt_.
        % - PriorTrunc   [double]   scalar, probability mass to be excluded, aka options_.prior_trunc
        % - Uniform      [logical]  scalar, produce uniform random deviates on the prior support.
        %
        % OUTPUTS
        % - o            [dprior]   scalar, prior object.
        %
        % REQUIREMENTS
        % None.
            if ~nargin % Empty object
                return
            end
            if isfield(bayestopt_, 'p1'), o.p1 = bayestopt_.p1; end
            if isfield(bayestopt_, 'p2'), o.p2 = bayestopt_.p2; end
            if isfield(bayestopt_, 'p3'), o.p3 = bayestopt_.p3; end
            if isfield(bayestopt_, 'p4'), o.p4 = bayestopt_.p4; end
            if isfield(bayestopt_, 'p5'), o.p5 = bayestopt_.p5; end
            if isfield(bayestopt_, 'p6'), o.p6 = bayestopt_.p6; end
            if isfield(bayestopt_, 'p7'), o.p7 = bayestopt_.p7; end
            if isfield(bayestopt_, 'p11'), o.p11 = bayestopt_.p11; end
            bounds = prior_bounds(bayestopt_, PriorTrunc);
            o.lb = bounds.lb;
            o.ub = bounds.ub;
            if nargin>2 && Uniform
                prior_shape = repmat(5, length(o.p6), 1);
            else
                prior_shape = bayestopt_.pshape;
            end
            o.idbeta = find(prior_shape==1);
            if ~isempty(o.idbeta)
                o.isbeta = true;
            end
            o.idgamma = find(prior_shape==2);
            if ~isempty(o.idgamma)
                o.isgamma = true;
            end
            o.idgaussian = find(prior_shape==3);
            if ~isempty(o.idgaussian)
                o.isgaussian = true;
            end
            o.idinvgamma1 = find(prior_shape==4);
            if ~isempty(o.idinvgamma1)
                o.isinvgamma1 = true;
            end
            o.iduniform = find(prior_shape==5);
            if ~isempty(o.iduniform)
                o.isuniform = true;
            end
            o.idinvgamma2 = find(prior_shape==6);
            if ~isempty(o.idinvgamma2)
                o.isinvgamma2 = true;
            end
            o.idweibull = find(prior_shape==8);
            if ~isempty(o.idweibull)
                o.isweibull = true;
            end
        end % dprior (constructor)

    end % methods

end % classdef --*-- Unit tests --*--

%@test:1
%$ % Fill global structures with required fields...
%$ prior_trunc = 1e-10;
%$ p0 = repmat([1; 2; 3; 4; 5; 6; 8], 2, 1);    % Prior shape
%$ p1 = .4*ones(14,1);                          % Prior mean
%$ p2 = .2*ones(14,1);                          % Prior std.
%$ p3 = NaN(14,1);
%$ p4 = NaN(14,1);
%$ p5 = NaN(14,1);
%$ p6 = NaN(14,1);
%$ p7 = NaN(14,1);
%$
%$ for i=1:14
%$    switch p0(i)
%$      case 1
%$        % Beta distribution
%$        p3(i) = 0;
%$        p4(i) = 1;
%$        [p6(i), p7(i)] = beta_specification(p1(i), p2(i)^2, p3(i), p4(i));
%$        p5(i) = compute_prior_mode([p6(i) p7(i)], 1);
%$      case 2
%$        % Gamma distribution
%$        p3(i) = 0;
%$        p4(i) = Inf;
%$        [p6(i), p7(i)] = gamma_specification(p1(i), p2(i)^2, p3(i), p4(i));
%$        p5(i) = compute_prior_mode([p6(i) p7(i)], 2);
%$      case 3
%$        % Normal distribution
%$        p3(i) = -Inf;
%$        p4(i) = Inf;
%$        p6(i) = p1(i);
%$        p7(i) = p2(i);
%$        p5(i) = p1(i);
%$      case 4
%$        % Inverse Gamma (type I) distribution
%$        p3(i) = 0;
%$        p4(i) = Inf;
%$        [p6(i), p7(i)] = inverse_gamma_specification(p1(i), p2(i)^2, p3(i), 1, false);
%$        p5(i) = compute_prior_mode([p6(i) p7(i)], 4);
%$      case 5
%$        % Uniform distribution
%$        [p1(i), p2(i), p6(i), p7(i)] = uniform_specification(p1(i), p2(i), p3(i), p4(i));
%$        p3(i) = p6(i);
%$        p4(i) = p7(i);
%$        p5(i) = compute_prior_mode([p6(i) p7(i)], 5);
%$      case 6
%$        % Inverse Gamma (type II) distribution
%$        p3(i) = 0;
%$        p4(i) = Inf;
%$        [p6(i), p7(i)] = inverse_gamma_specification(p1(i), p2(i)^2, p3(i), 2, false);
%$        p5(i) = compute_prior_mode([p6(i) p7(i)], 6);
%$      case 8
%$        % Weibull distribution
%$        p3(i) = 0;
%$        p4(i) = Inf;
%$        [p6(i), p7(i)] = weibull_specification(p1(i), p2(i)^2, p3(i));
%$        p5(i) = compute_prior_mode([p6(i) p7(i)], 8);
%$      otherwise
%$        error('This density is not implemented!')
%$    end
%$ end
%$
%$ bayestopt_.pshape = p0;
%$ bayestopt_.p1 = p1;
%$ bayestopt_.p2 = p2;
%$ bayestopt_.p3 = p3;
%$ bayestopt_.p4 = p4;
%$ bayestopt_.p5 = p5;
%$ bayestopt_.p6 = p6;
%$ bayestopt_.p7 = p7;
%$
%$ ndraws = 1e5;
%$ m0 = bayestopt_.p1; %zeros(14,1);
%$ v0 = diag(bayestopt_.p2.^2); %zeros(14);
%$
%$ % Call the tested routine
%$ try
%$    % Instantiate dprior object
%$    o = dprior(bayestopt_, prior_trunc, false);
%$    t(1) = true;
%$ catch
%$     t(1) = false;
%$ end
%$
%$ T = all(t);
%@eof:1

%@test:2
%$ try
%$    % Instantiate dprior object
%$    o = dprior();
%$    t(1) = true;
%$ catch
%$     t(1) = false;
%$ end
%$
%$ T = all(t);
%@eof:2
