function mcheck = mode_check(fun,xparam,hessian_mat,options_,M_,estim_params_,bayestopt_,bounds,isMinimum, varargin)
% function mcheck = mode_check(fun,xparam,hessian_mat,options_,M_,estim_params_,bayestopt_,bounds,isMinimum, varargin)
% -------------------------------------------------------------------------
% Checks the estimated ML or Posterior mode/minimum by plotting sections of
% the likelihood/posterior kernel. Each plot shows the variation of the
% function implied by the variations of a single parameter ( ceteris paribus)
% -------------------------------------------------------------------------
% INPUTS
% - fun:            [func_handle]  objective function
% - xparam:         [vector]       estimated mode/minimum
% - hessian_mat:    [matrix]       hessian of the objective function at the estimated mode/minimum
% - options_:       [structure]    Dynare options structure
% - M_:             [structure]    Dynare model structure
% - estim_params_:  [structure]    Dynare estimated parameters structure
% - bayestopt_:     [structure]    information on the priors
% - bounds:         [structure]    information on the bounds
% - isMinimum:      [boolean]      true if xparam is a minimum, false if it is a mode
% - varargin:       [cell]         additional arguments to be passed to fun
% -------------------------------------------------------------------------
% OUTPUTS
% - mcheck: [structure]     structure containing the data for the check plots
% - Saves the plots in the graphs folder and the data in a mat file
% -------------------------------------------------------------------------
% This function is called by
% - dynare_estimation_1
% - mom.run
% -------------------------------------------------------------------------

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

tolBounds = 1e-8;

fval = feval(fun,xparam,varargin{:});

if ~isempty(hessian_mat)
    [ s_min, k ] = min(diag(hessian_mat));
    if isMinimum
        fprintf('\nMINIMUM CHECK\n\nFval obtained by the optimization routine: %f\n', fval)
    else
        fprintf('\nMODE CHECK\n\nFval obtained by the optimization routine: %f\n', fval)        
    end
    if s_min<eps
        fprintf('Most negative variance %f for parameter %d (%s = %f)\n', s_min, k , bayestopt_.name{k}, xparam(k));
    end
end

[nbplt,nr,nc,~,~,nstar] = pltorg(length(xparam));

graphsFolder = CheckPath('graphs',M_.dname);
latexFolder = CheckPath('latex',M_.dname);

if options_.TeX && any(strcmp('eps',cellstr(options_.graph_format)))
    fidTeX = fopen([latexFolder filesep M_.fname '_CheckPlots.tex'],'w');
    fprintf(fidTeX,'%% TeX eps-loader file generated by mode_check.m (Dynare).\n');
    fprintf(fidTeX,['%% ' datestr(now,0) '\n']);
    fprintf(fidTeX,' \n');
end

ll = options_.mode_check.neighbourhood_size;
if isinf(ll)
    options_.mode_check.symmetric_plots = false;
end

if isMinimum
    mcheck = struct('cross',struct(),'emin',struct());
else
    mcheck = struct('cross',struct(),'emode',struct());
end

for plt = 1:nbplt
    if options_.TeX
        NAMES = [];
        TeXNAMES = [];
    end
    if isMinimum
        hh_fig = dyn_figure(options_.nodisplay,'Name','Minimum check plots');
    else
        hh_fig = dyn_figure(options_.nodisplay,'Name','Mode check plots');
    end
    for k = 1:min(nstar,length(xparam)-(plt-1)*nstar)
        subplot(nr,nc,k)
        kk = (plt-1)*nstar+k;
        [name,texname] = get_the_name(kk,options_.TeX,M_,estim_params_,options_.varobs);
        xx = xparam;
        if xparam(kk)~=0 && ~isinf(bounds.lb(kk)) && ~isinf(bounds.ub(kk))
            l1 = max(bounds.lb(kk),(1-sign(xparam(kk))*ll)*xparam(kk)); m1 = 0; % lower bound
            l2 = min(bounds.ub(kk),(1+sign(xparam(kk))*ll)*xparam(kk));         % upper bound
        else
            % size info for 0 parameter is missing, use prior standard deviation
            upper_bound = bounds.lb(kk);
            if isinf(upper_bound)
                upper_bound = -1e-6*options_.huge_number;
            end
            lower_bound = bounds.ub(kk);
            if isinf(lower_bound)
                lower_bound = -1e-6*options_.huge_number;
            end
            l1 = max(lower_bound,-bayestopt_.p2(kk)); m1 = 0; % lower bound
            l2 = min(upper_bound,bayestopt_.p2(kk));          % upper bound
        end
        binding_lower_bound = 0;
        binding_upper_bound = 0;
        if abs(xparam(kk)-bounds.lb(kk))<tolBounds
            binding_lower_bound = 1;
            bound_value = bounds.lb(kk);
        elseif abs(xparam(kk)-bounds.ub(kk))<tolBounds
            binding_upper_bound = 1;
            bound_value = bounds.ub(kk);
        end
        if options_.mode_check.symmetric_plots && ~binding_lower_bound && ~binding_upper_bound
            if l2<(1+ll)*xparam(kk) % test whether upper bound is too small due to prior binding
                l1 = xparam(kk) - (l2-xparam(kk)); % adjust lower bound to become closer
                m1 = 1;
            end
            if ~m1 && (l1>(1-ll)*xparam(kk)) && (xparam(kk)+(xparam(kk)-l1)<bounds.ub(kk)) % if lower bound was truncated and using difference from lower bound does not violate upper bound
                l2 = xparam(kk) + (xparam(kk)-l1); % set upper bound to same distance as lower bound
            end
        end
        z1 = l1:((xparam(kk)-l1)/(options_.mode_check.number_of_points/2)):xparam(kk);
        z2 = xparam(kk):((l2-xparam(kk))/(options_.mode_check.number_of_points/2)):l2;
        z  = union(z1,z2);
        if ~options_.mode_check.nolik
            y = zeros(length(z),2);
            if isfield(options_,'mom') && ( (strcmp(options_.mom.mom_method,'GMM') || strcmp(options_.mom.mom_method,'SMM')) && options_.mom.penalized_estimator )
                dy = (xx-bayestopt_.p1)'/diag(bayestopt_.p2.^2)*(xx-bayestopt_.p1);
            else
                dy = priordens(xx,bayestopt_.pshape,bayestopt_.p6,bayestopt_.p7,bayestopt_.p3,bayestopt_.p4);
            end
        end
        for i = 1:length(z)
            xx(kk) = z(i);
            [fval, info, exit_flag] = feval(fun,xx, varargin{:});                
            if exit_flag
                y(i,1) = fval;
            else
                y(i,1) = NaN;
                if options_.debug
                    fprintf('mode_check:: could not solve model for parameter %s at value %4.3f, error code: %u (%s)\n',name,z(i),info(1),get_error_message(info, options_));
                end
            end
            if ~options_.mode_check.nolik
                if isfield(options_,'mom') && ( (strcmp(options_.mom.mom_method,'GMM') || strcmp(options_.mom.mom_method,'SMM')) && options_.mom.penalized_estimator )
                    lnprior = (xx-bayestopt_.p1)'/diag(bayestopt_.p2.^2)*(xx-bayestopt_.p1);
                else
                    lnprior = priordens(xx,bayestopt_.pshape,bayestopt_.p6,bayestopt_.p7,bayestopt_.p3,bayestopt_.p4);
                end
                y(i,2)  = (y(i,1)+lnprior-dy);
            end
        end
        if isMinimum
            mcheck.cross = setfield(mcheck.cross, name, [transpose(z), y]); % keep y
            mcheck.emin = setfield(mcheck.emin, name, xparam(kk));          % store as min
            fighandle = plot(z,y);
        else
            mcheck.cross = setfield(mcheck.cross, name, [transpose(z), -y]); % multiply y by -1
            mcheck.emode = setfield(mcheck.emode, name, xparam(kk));         % store as mode
            fighandle = plot(z,-y);
        end        
        hold on
        yl = get(gca,'ylim');
        plot( [xparam(kk) xparam(kk)], yl, 'c', 'LineWidth', 1);
        NaN_index = find(isnan(y(:,1)));
        zNaN = z(NaN_index);
        yNaN = yl(1)*ones(size(NaN_index));
        plot(zNaN,yNaN,'o','MarkerEdgeColor','r','MarkerFaceColor','r','MarkerSize',6);
        if options_.TeX
            title(texname,'interpreter','latex')
        else
            title(name,'interpreter','none')
        end
        axis tight
        if binding_lower_bound || binding_upper_bound
            xl = get(gca,'xlim');
            plot( [bound_value bound_value], yl, 'r--', 'LineWidth', 1);
            xlim([xl(1)-0.5*binding_lower_bound*(xl(2)-xl(1)) xl(2)+0.5*binding_upper_bound*(xl(2)-xl(1))]);
        end
        hold off
        drawnow
    end
    if ~options_.mode_check.nolik
        if isoctave
            axes('outerposition',[0.3 0.93 0.42 0.07],'box','on');
        else
            axes('position',[0.3 0.01 0.42 0.05],'box','on');
        end
        line_color=get(fighandle,'color');
        plot([0.48 0.68],[0.5 0.5],'color',line_color{2});
        hold on;
        plot([0.04 0.24],[0.5 0.5],'color',line_color{1});
        set(gca,'xlim',[0 1],'ylim',[0 1],'xtick',[],'ytick',[]);
        if isMinimum
            text(0.25,0.5,'log-post');
            text(0.69,0.5,'log-dist kernel');
        else
            text(0.25,0.5,'log-post');
            text(0.69,0.5,'log-lik kernel');
        end
    end
    dyn_saveas(hh_fig,[graphsFolder filesep M_.fname '_CheckPlots' int2str(plt) ],options_.nodisplay,options_.graph_format);
    if options_.TeX && any(strcmp('eps',cellstr(options_.graph_format)))
        % TeX eps loader file
        fprintf(fidTeX,'\\begin{figure}[H]\n');
        fprintf(fidTeX,'\\centering \n');
        fprintf(fidTeX,'\\includegraphics[width=%2.2f\\textwidth]{%s_CheckPlots%s}\n',options_.figures.textwidth*min(k/nc,1),[graphsFolder '/' M_.fname],int2str(plt)); % don't use filesep as it will create issues with LaTeX on Windows
        fprintf(fidTeX,'\\caption{Check plots.}');
        fprintf(fidTeX,'\\label{Fig:CheckPlots:%s}\n',int2str(plt));
        fprintf(fidTeX,'\\end{figure}\n');
        fprintf(fidTeX,' \n');
    end
end
if options_.TeX && any(strcmp('eps',cellstr(options_.graph_format)))
    fclose(fidTeX);
end

save([graphsFolder filesep M_.fname '_check_plot_data.mat'],'mcheck');