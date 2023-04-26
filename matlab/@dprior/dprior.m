classdef dprior

    properties
        p6 = [];                         % Prior first hyperparameter.
        p7 = [];                         % Prior second hyperparameter.
        p3 = [];                         % Lower bound of the prior support.
        p4 = [];                         % Upper bound of the prior support.
        lb = [];                         % Truncated prior lower bound.
        ub = [];                         % Truncated prior upper bound.
        uniform_index = [];              % Index for the uniform priors.
        gaussian_index = [];             % Index for the gaussian priors.
        gamma_index = [];                % Index for the gamma priors.
        beta_index = [];                 % Index for the beta priors.
        inverse_gamma_1_index = [];      % Index for the inverse gamma type 1 priors.
        inverse_gamma_2_index = [];      % Index for the inverse gamma type 2 priors.
        weibull_index = [];              % Index for the weibull priors.
        uniform_draws = false;
        gaussian_draws = false;
        gamma_draws = false;
        beta_draws = false;
        inverse_gamma_1_draws = false;
        inverse_gamma_2_draws = false;
        weibull_draws = false;
    end

    methods

        function o = dprior(BayesInfo, PriorTrunc, Uniform)
        % Class constructor.
        %
        % INPUTS
        % - BayesInfo    [struct]   Informations about the prior distribution, aka bayestopt_.
        % - PriorTrunc   [double]   scalar, probability mass to be excluded, aka options_.prior_trunc
        % - Uniform      [logical]  scalar, produce uniform random deviates on the prior support.
        %
        % OUTPUTS
        % - o            [dprior]   scalar, prior object.
        %
        % REQUIREMENTS
        % None.
            o.p6 = BayesInfo.p6;
            o.p7 = BayesInfo.p7;
            o.p3 = BayesInfo.p3;
            o.p4 = BayesInfo.p4;
            bounds = prior_bounds(BayesInfo, PriorTrunc);
            o.lb = bounds.lb;
            o.ub = bounds.ub;
            if nargin>2 && Uniform
                prior_shape = repmat(5, length(o.p6), 1);
            else
                prior_shape = BayesInfo.pshape;
            end
            o.beta_index = find(prior_shape==1);
            if ~isempty(o.beta_index)
                o.beta_draws = true;
            end
            o.gamma_index = find(prior_shape==2);
            if ~isempty(o.gamma_index)
                o.gamma_draws = true;
            end
            o.gaussian_index = find(prior_shape==3);
            if ~isempty(o.gaussian_index)
                o.gaussian_draws = true;
            end
            o.inverse_gamma_1_index = find(prior_shape==4);
            if ~isempty(o.inverse_gamma_1_index)
                o.inverse_gamma_1_draws = true;
            end
            o.uniform_index = find(prior_shape==5);
            if ~isempty(o.uniform_index)
                o.uniform_draws = true;
            end
            o.inverse_gamma_2_index = find(prior_shape==6);
            if ~isempty(o.inverse_gamma_2_index)
                o.inverse_gamma_2_draws = true;
            end
            o.weibull_index = find(prior_shape==8);
            if ~isempty(o.weibull_index)
                o.weibull_draws = true;
            end
        end

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
            p = NaN(rows(o.lb), 1);
            if o.uniform_draws
                p(o.uniform_index) = rand(length(o.uniform_index),1).*(o.p4(o.uniform_index)-o.p3(o.uniform_index)) + o.p3(o.uniform_index);
                out_of_bound = find( (p(o.uniform_index)>o.ub(o.uniform_index)) | (p(o.uniform_index)<o.lb(o.uniform_index)));
                while ~isempty(out_of_bound)
                    p(o.uniform_index) = rand(length(o.uniform_index), 1).*(o.p4(o.uniform_index)-o.p3(o.uniform_index)) + o.p3(o.uniform_index);
                    out_of_bound = find( (p(o.uniform_index)>o.ub(o.uniform_index)) | (p(o.uniform_index)<o.lb(o.uniform_index)));
                end
            end
            if o.gaussian_draws
                p(o.gaussian_index) = randn(length(o.gaussian_index), 1).*o.p7(o.gaussian_index) + o.p6(o.gaussian_index);
                out_of_bound = find( (p(o.gaussian_index)>o.ub(o.gaussian_index)) | (p(o.gaussian_index)<o.lb(o.gaussian_index)));
                while ~isempty(out_of_bound)
                    p(o.gaussian_index(out_of_bound)) = randn(length(o.gaussian_index(out_of_bound)), 1).*o.p7(o.gaussian_index(out_of_bound)) + o.p6(o.gaussian_index(out_of_bound));
                    out_of_bound = find( (p(o.gaussian_index)>o.ub(o.gaussian_index)) | (p(o.gaussian_index)<o.lb(o.gaussian_index)));
                end
            end
            if o.gamma_draws
                p(o.gamma_index) = gamrnd(o.p6(o.gamma_index), o.p7(o.gamma_index))+o.p3(o.gamma_index);
                out_of_bound = find( (p(o.gamma_index)>o.ub(o.gamma_index)) | (p(o.gamma_index)<o.lb(o.gamma_index)));
                while ~isempty(out_of_bound)
                    p(o.gamma_index(out_of_bound)) = gamrnd(o.p6(o.gamma_index(out_of_bound)), o.p7(o.gamma_index(out_of_bound)))+o.p3(o.gamma_index(out_of_bound));
                    out_of_bound = find( (p(o.gamma_index)>o.ub(o.gamma_index)) | (p(o.gamma_index)<o.lb(o.gamma_index)));
                end
            end
            if o.beta_draws
                p(o.beta_index) = (o.p4(o.beta_index)-o.p3(o.beta_index)).*betarnd(o.p6(o.beta_index), o.p7(o.beta_index))+o.p3(o.beta_index);
                out_of_bound = find( (p(o.beta_index)>o.ub(o.beta_index)) | (p(o.beta_index)<o.lb(o.beta_index)));
                while ~isempty(out_of_bound)
                    p(o.beta_index(out_of_bound)) = (o.p4(o.beta_index(out_of_bound))-o.p3(o.beta_index(out_of_bound))).*betarnd(o.p6(o.beta_index(out_of_bound)), o.p7(o.beta_index(out_of_bound)))+o.p3(o.beta_index(out_of_bound));
                    out_of_bound = find( (p(o.beta_index)>o.ub(o.beta_index)) | (p(o.beta_index)<o.lb(o.beta_index)));
                end
            end
            if o.inverse_gamma_1_draws
                p(o.inverse_gamma_1_index) = ...
                    sqrt(1./gamrnd(o.p7(o.inverse_gamma_1_index)/2, 2./o.p6(o.inverse_gamma_1_index)))+o.p3(o.inverse_gamma_1_index);
                out_of_bound = find( (p(o.inverse_gamma_1_index)>o.ub(o.inverse_gamma_1_index)) | (p(o.inverse_gamma_1_index)<o.lb(o.inverse_gamma_1_index)));
                while ~isempty(out_of_bound)
                    p(o.inverse_gamma_1_index(out_of_bound)) = ...
                        sqrt(1./gamrnd(o.p7(o.inverse_gamma_1_index(out_of_bound))/2, 2./o.p6(o.inverse_gamma_1_index(out_of_bound))))+o.p3(o.inverse_gamma_1_index(out_of_bound));
                    out_of_bound = find( (p(o.inverse_gamma_1_index)>o.ub(inverse_gamma_1_index)) | (p(o.inverse_gamma_1_index)<o.lb(o.inverse_gamma_1_index)));
                end
            end
            if o.inverse_gamma_2_draws
                p(o.inverse_gamma_2_index) = ...
                    1./gamrnd(o.p7(o.inverse_gamma_2_index)/2, 2./o.p6(o.inverse_gamma_2_index))+o.p3(o.inverse_gamma_2_index);
                out_of_bound = find( (p(o.inverse_gamma_2_index)>o.ub(o.inverse_gamma_2_index)) | (p(o.inverse_gamma_2_index)<o.lb(o.inverse_gamma_2_index)));
                while ~isempty(out_of_bound)
                    p(o.inverse_gamma_2_index(out_of_bound)) = ...
                        1./gamrnd(o.p7(o.inverse_gamma_2_index(out_of_bound))/2, 2./o.p6(o.inverse_gamma_2_index(out_of_bound)))+o.p3(o.inverse_gamma_2_index(out_of_bound));
                    out_of_bound = find( (p(o.inverse_gamma_2_index)>o.ub(o.inverse_gamma_2_index)) | (p(o.inverse_gamma_2_index)<o.lb(o.inverse_gamma_2_index)));
                end
            end
            if o.weibull_draws
                p(o.weibull_index) = wblrnd(o.p7(o.weibull_index), o.p6(o.weibull_index)) + o.p3(o.weibull_index);
                out_of_bound = find( (p(o.weibull_index)>o.ub(o.weibull_index)) | (p(o.weibull_index)<o.lb(o.weibull_index)));
                while ~isempty(out_of_bound)
                    p(o.weibull_index(out_of_bound)) = wblrnd(o.p7(o.weibull_index(out_of_bound)), o.p6(o.weibull_index(out_of_bound)))+o.p3(o.weibull_index(out_of_bound));
                    out_of_bound = find( (p(o.weibull_index)>o.ub(o.weibull_index)) | (p(o.weibull_index)<o.lb(o.weibull_index)));
                end
            end
        end

        function P = draws(o, n)
        % Return n independent random draws from the prior distribution.
        %
        % INPUTS
        % - o    [dprior]
        %
        % OUTPUTS
        % - P    [double]   m×n matrix, random draw from the prior distribution.
        %
        % REMARKS
        % If the Parallel Computing Toolbox is available, the main loop is run in parallel.
        %
        % EXAMPLE
        %
        % >> Prior = dprior(bayestopt_, options_.prior_trunc);
        % >> Prior.draws(1e6)
            P = NaN(rows(o.lb), 1);
            parfor i=1:n
                P(:,i) = draw(o);
            end
        end

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
%$ BayesInfo.pshape = p0;
%$ BayesInfo.p1 = p1;
%$ BayesInfo.p2 = p2;
%$ BayesInfo.p3 = p3;
%$ BayesInfo.p4 = p4;
%$ BayesInfo.p5 = p5;
%$ BayesInfo.p6 = p6;
%$ BayesInfo.p7 = p7;
%$
%$ ndraws = 1e5;
%$ m0 = BayesInfo.p1; %zeros(14,1);
%$ v0 = diag(BayesInfo.p2.^2); %zeros(14);
%$
%$ % Call the tested routine
%$ try
%$    % Instantiate dprior object
%$    o = dprior(BayesInfo, prior_trunc, false);
%$    % Do simulations in a loop and estimate recursively the mean and the variance.
%$    for i = 1:ndraws
%$         draw = o.draw();
%$         m1 = m0 + (draw-m0)/i;
%$         m2 = m1*m1';
%$         v0 = v0 + ((draw*draw'-m2-v0) + (i-1)*(m0*m0'-m2'))/i;
%$         m0 = m1;
%$    end
%$    t(1) = true;
%$ catch
%$     t(1) = false;
%$ end
%$
%$ if t(1)
%$     t(2) = all(abs(m0-BayesInfo.p1)<3e-3);
%$     t(3) = all(all(abs(v0-diag(BayesInfo.p2.^2))<5e-3));
%$ end
%$ T = all(t);
%@eof:1

%@test:2
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
%$ BayesInfo.pshape = p0;
%$ BayesInfo.p1 = p1;
%$ BayesInfo.p2 = p2;
%$ BayesInfo.p3 = p3;
%$ BayesInfo.p4 = p4;
%$ BayesInfo.p5 = p5;
%$ BayesInfo.p6 = p6;
%$ BayesInfo.p7 = p7;
%$
%$ ndraws = 1e5;
%$
%$ % Call the tested routine
%$ try
%$    % Instantiate dprior object.
%$    o = dprior(BayesInfo, prior_trunc, false);
%$    X = o.draws(ndraws);
%$    m = mean(X, 2);
%$    v = var(X, 0, 2);
%$    t(1) = true;
%$ catch
%$     t(1) = false;
%$ end
%$
%$ if t(1)
%$     t(2) = all(abs(m-BayesInfo.p1)<3e-3);
%$     t(3) = all(all(abs(v-BayesInfo.p2.^2)<5e-3));
%$ end
%$ T = all(t);
%@eof:2
