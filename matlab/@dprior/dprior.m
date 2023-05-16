classdef dprior < handle

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

        function p = subsref(o, S)
            switch S(1).type
              case '.'
                if ismember(S(1).subs, {'p1','p2','p3','p4','p5','p6','p7','lb','ub'})
                    p = builtin('subsref', o, S(1));
                elseif ismember(S(1).subs, {'draw'})
                    p = feval(S(1).subs, o);
                elseif ismember(S(1).subs, {'draws', 'density', 'densities', 'moments'})
                    p = feval(S(1).subs, o , S(2).subs{:});
                elseif ismember(S(1).subs, {'mean', 'median', 'variance', 'mode'})
                    if (length(S)==2 && isempty(S(2).subs)) || length(S)==1
                        p = feval(S(1).subs, o);
                    else
                        p = feval(S(1).subs, o , S(2).subs{:});
                    end
                else
                    error('dprior::subsref: unknown method (%s).', S(1).subs)
                end
              otherwise
                error('dprior::subsref: %s indexing not implemented.', S(1).type)
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
        end % draw

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
        end % draws

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
        end % density

        function lpd = densities(o, X)
        % Evaluate the logged prior densities at X (for each column).
        %
        % INPUTS
        % - o       [dprior]
        % - X       [double]   m×n matrix, n points where the prior density is evaluated.
        %
        % OUTPUTS
        % - lpd     [double]   1×n, values of the logged prior density at X.
            n = columns(X);
            lpd = NaN(1, n);
            parfor i=1:n
                lpd(i) = density(o, X(:,i));
            end
        end % densities

        function o = moments(o, name)
        % Compute the prior moments.
        %
        % INPUTS
        % - o       [dprior]
        %
        % OUTPUTS
        % - o       [dprior]
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
        end

        function m = mean(o, resetmoments)
        % Return the prior mean.
        %
        % INPUTS
        % - o               [dprior]
        % - resetmoments    [logical]     Force the computation of the prior mean
        %
        % OUTPUTS
        % - m               [double]      n×1 vector, prior mean
            if nargin<2, resetmoments = false; end
            if any(isnan(o.p1)), resetmoments = true; end
            if resetmoments, o.p1 = NaN(size(o.p1)); o.moments('mean');
            end
                m = o.p1;
            end

            function m = variance(o, resetmoments)
            % Return the prior variance.
            %
            % INPUTS
            % - o               [dprior]
            % - resetmoments    [logical]     Force the computation of the prior variance
            %
            % OUTPUTS
            % - m               [double]      n×1 vector, prior variance
                if nargin<2, resetmoments = false; end
                if any(isnan(o.p2)), resetmoments = true; end
                if resetmoments, o.p2 = NaN(size(o.p2)); o.moments('std'); end
                m = o.p2.^2;
            end

            function m = median(o, resetmoments)
            % Return the prior median.
            %
            % INPUTS
            % - o               [dprior]
            % - resetmoments    [logical]     Force the computation of the prior median
            %
            % OUTPUTS
            % - m               [double]      n×1 vector, prior median
                if nargin<2, resetmoments = false; end
                if any(isnan(o.p11)), resetmoments = true; end
                if resetmoments, o.p11 = NaN(size(o.p11)); o.moments('median'); end
                m = o.p11;
            end

            function m = mode(o, resetmoments)
            % Return the prior mode.
            %
            % INPUTS
            % - o               [dprior]
            % - resetmoments    [logical]     Force the computation of the prior mode
            %
            % OUTPUTS
            % - m               [double]      n×1 vector, prior mode
                if nargin<2, resetmoments = false; end
                if any(isnan(o.p5)), resetmoments = true; end
                if resetmoments, o.p5 = NaN(size(o.p5)); o.moments('mode'); end
                m = o.p5;
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
%$     t(2) = all(abs(m0-bayestopt_.p1)<3e-3);
%$     t(3) = all(all(abs(v0-diag(bayestopt_.p2.^2))<5e-3));
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
%$
%$ % Call the tested routine
%$ try
%$    % Instantiate dprior object.
%$    o = dprior(bayestopt_, prior_trunc, false);
%$    X = o.draws(ndraws);
%$    m = mean(X, 2);
%$    v = var(X, 0, 2);
%$    t(1) = true;
%$ catch
%$     t(1) = false;
%$ end
%$
%$ if t(1)
%$     t(2) = all(abs(m-bayestopt_.p1)<3e-3);
%$     t(3) = all(all(abs(v-bayestopt_.p2.^2)<5e-3));
%$ end
%$ T = all(t);
%@eof:2

%@test:3
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
%$     switch p0(i)
%$       case 1
%$         % Beta distribution
%$         p3(i) = 0;
%$         p4(i) = 1;
%$         [p6(i), p7(i)] = beta_specification(p1(i), p2(i)^2, p3(i), p4(i));
%$         p5(i) = compute_prior_mode([p6(i) p7(i)], 1);
%$       case 2
%$         % Gamma distribution
%$         p3(i) = 0;
%$         p4(i) = Inf;
%$         [p6(i), p7(i)] = gamma_specification(p1(i), p2(i)^2, p3(i), p4(i));
%$         p5(i) = compute_prior_mode([p6(i) p7(i)], 2);
%$       case 3
%$         % Normal distribution
%$         p3(i) = -Inf;
%$         p4(i) = Inf;
%$         p6(i) = p1(i);
%$         p7(i) = p2(i);
%$         p5(i) = p1(i);
%$       case 4
%$         % Inverse Gamma (type I) distribution
%$         p3(i) = 0;
%$         p4(i) = Inf;
%$         [p6(i), p7(i)] = inverse_gamma_specification(p1(i), p2(i)^2, p3(i), 1, false);
%$         p5(i) = compute_prior_mode([p6(i) p7(i)], 4);
%$       case 5
%$         % Uniform distribution
%$         [p1(i), p2(i), p6(i), p7(i)] = uniform_specification(p1(i), p2(i), p3(i), p4(i));
%$         p3(i) = p6(i);
%$         p4(i) = p7(i);
%$         p5(i) = compute_prior_mode([p6(i) p7(i)], 5);
%$       case 6
%$         % Inverse Gamma (type II) distribution
%$         p3(i) = 0;
%$         p4(i) = Inf;
%$         [p6(i), p7(i)] = inverse_gamma_specification(p1(i), p2(i)^2, p3(i), 2, false);
%$         p5(i) = compute_prior_mode([p6(i) p7(i)], 6);
%$       case 8
%$         % Weibull distribution
%$         p3(i) = 0;
%$         p4(i) = Inf;
%$         [p6(i), p7(i)] = weibull_specification(p1(i), p2(i)^2, p3(i));
%$         p5(i) = compute_prior_mode([p6(i) p7(i)], 8);
%$       otherwise
%$         error('This density is not implemented!')
%$     end
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
%$ % Call the tested routine
%$ try
%$     Prior = dprior(bayestopt_, prior_trunc, false);
%$
%$     % Compute density at the prior mode
%$     lpdstar = Prior.density(p5);
%$
%$     % Draw random deviates in a loop and evaluate the density.
%$     LPD = NaN(10000,1);
%$     parfor i = 1:10000
%$         x = Prior.draw();
%$         LPD(i) = Prior.density(x);
%$     end
%$     t(1) = true;
%$ catch
%$     t(1) = false;
%$ end
%$
%$ if t(1)
%$     t(2) = all(LPD<=lpdstar);
%$ end
%$ T = all(t);
%@eof:3

%@test:4
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
%$     switch p0(i)
%$       case 1
%$         % Beta distribution
%$         p3(i) = 0;
%$         p4(i) = 1;
%$         [p6(i), p7(i)] = beta_specification(p1(i), p2(i)^2, p3(i), p4(i));
%$         p5(i) = compute_prior_mode([p6(i) p7(i)], 1);
%$       case 2
%$         % Gamma distribution
%$         p3(i) = 0;
%$         p4(i) = Inf;
%$         [p6(i), p7(i)] = gamma_specification(p1(i), p2(i)^2, p3(i), p4(i));
%$         p5(i) = compute_prior_mode([p6(i) p7(i)], 2);
%$       case 3
%$         % Normal distribution
%$         p3(i) = -Inf;
%$         p4(i) = Inf;
%$         p6(i) = p1(i);
%$         p7(i) = p2(i);
%$         p5(i) = p1(i);
%$       case 4
%$         % Inverse Gamma (type I) distribution
%$         p3(i) = 0;
%$         p4(i) = Inf;
%$         [p6(i), p7(i)] = inverse_gamma_specification(p1(i), p2(i)^2, p3(i), 1, false);
%$         p5(i) = compute_prior_mode([p6(i) p7(i)], 4);
%$       case 5
%$         % Uniform distribution
%$         [p1(i), p2(i), p6(i), p7(i)] = uniform_specification(p1(i), p2(i), p3(i), p4(i));
%$         p3(i) = p6(i);
%$         p4(i) = p7(i);
%$         p5(i) = compute_prior_mode([p6(i) p7(i)], 5);
%$       case 6
%$         % Inverse Gamma (type II) distribution
%$         p3(i) = 0;
%$         p4(i) = Inf;
%$         [p6(i), p7(i)] = inverse_gamma_specification(p1(i), p2(i)^2, p3(i), 2, false);
%$         p5(i) = compute_prior_mode([p6(i) p7(i)], 6);
%$       case 8
%$         % Weibull distribution
%$         p3(i) = 0;
%$         p4(i) = Inf;
%$         [p6(i), p7(i)] = weibull_specification(p1(i), p2(i)^2, p3(i));
%$         p5(i) = compute_prior_mode([p6(i) p7(i)], 8);
%$       otherwise
%$         error('This density is not implemented!')
%$     end
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
%$ % Call the tested routine
%$ try
%$     Prior = dprior(bayestopt_, prior_trunc, false);
%$     mu = NaN(14,1);
%$     std = NaN(14,1);
%$
%$     for i=1:14
%$         % Define conditional density (it's also a marginal since priors are orthogonal)
%$         f = @(x) exp(Prior.densities(substitute(p5, i, x)));
%$         % TODO: Check the version of Octave we use (integral is available as a compatibility wrapper in latest Octave version)
%$         m = integral(f, p3(i), p4(i));
%$         density = @(x) f(x)/m; % rescaling is required since the probability mass depends on the conditioning.
%$         % Compute the conditional expectation
%$         mu(i) = integral(@(x) x.*density(x), p3(i), p4(i));
%$         std(i) = sqrt(integral(@(x) ((x-mu(i)).^2).*density(x), p3(i), p4(i)));
%$     end
%$
%$     t(1) = true;
%$ catch
%$     t(1) = false;
%$ end
%$
%$ if t(1)
%$     t(2) = all(abs(mu-.4)<1e-6);
%$     t(3) = all(abs(std-.2)<1e-6);
%$ end
%$ T = all(t);
%@eof:4

%@test:5
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
%$     switch p0(i)
%$       case 1
%$         % Beta distribution
%$         p3(i) = 0;
%$         p4(i) = 1;
%$         [p6(i), p7(i)] = beta_specification(p1(i), p2(i)^2, p3(i), p4(i));
%$         p5(i) = compute_prior_mode([p6(i) p7(i)], 1);
%$       case 2
%$         % Gamma distribution
%$         p3(i) = 0;
%$         p4(i) = Inf;
%$         [p6(i), p7(i)] = gamma_specification(p1(i), p2(i)^2, p3(i), p4(i));
%$         p5(i) = compute_prior_mode([p6(i) p7(i)], 2);
%$       case 3
%$         % Normal distribution
%$         p3(i) = -Inf;
%$         p4(i) = Inf;
%$         p6(i) = p1(i);
%$         p7(i) = p2(i);
%$         p5(i) = p1(i);
%$       case 4
%$         % Inverse Gamma (type I) distribution
%$         p3(i) = 0;
%$         p4(i) = Inf;
%$         [p6(i), p7(i)] = inverse_gamma_specification(p1(i), p2(i)^2, p3(i), 1, false);
%$         p5(i) = compute_prior_mode([p6(i) p7(i)], 4);
%$       case 5
%$         % Uniform distribution
%$         [p1(i), p2(i), p6(i), p7(i)] = uniform_specification(p1(i), p2(i), p3(i), p4(i));
%$         p3(i) = p6(i);
%$         p4(i) = p7(i);
%$         p5(i) = compute_prior_mode([p6(i) p7(i)], 5);
%$       case 6
%$         % Inverse Gamma (type II) distribution
%$         p3(i) = 0;
%$         p4(i) = Inf;
%$         [p6(i), p7(i)] = inverse_gamma_specification(p1(i), p2(i)^2, p3(i), 2, false);
%$         p5(i) = compute_prior_mode([p6(i) p7(i)], 6);
%$       case 8
%$         % Weibull distribution
%$         p3(i) = 0;
%$         p4(i) = Inf;
%$         [p6(i), p7(i)] = weibull_specification(p1(i), p2(i)^2, p3(i));
%$         p5(i) = compute_prior_mode([p6(i) p7(i)], 8);
%$       otherwise
%$         error('This density is not implemented!')
%$     end
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
%$ % Call the tested routine
%$ try
%$     Prior = dprior(BayesInfo, prior_trunc, false);
%$     t(1) = true;
%$ catch
%$     t(1) = false;
%$ end
%$
%$ if t(1)
%$     t(2) = all(Prior.mean()==.4);
%$     t(3) = all(ismembertol(Prior.mean(true),.4));
%$     t(4) = all(ismembertol(Prior.variance(),.04));
%$     t(5) = all(ismembertol(Prior.variance(true),.04));
%$ end
%$ T = all(t);
%@eof:5
