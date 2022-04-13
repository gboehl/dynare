function [info, oo_, options_, M_] = stoch_simul(M_, options_, oo_, var_list)

% Copyright © 2001-2021 Dynare Team
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

% Test if the order of approximation is nonzero (the preprocessor tests if order is non negative).
if isequal(options_.order,0)
    error('stoch_simul:: The order of the Taylor approximation cannot be 0!')
end

% The TeX option crashes MATLAB R2014a run with "-nodisplay" option
% (as is done from the testsuite).
% Since we can’t directly test whether "-nodisplay" has been passed,
% we test for the "TOP_TEST_DIR" environment variable, which is set
% by the testsuite.
% Note that it was not tested whether the crash happens with more
% recent MATLAB versions, so when OLD_MATLAB_VERSION is increased,
% one should make a test before removing this workaround.
if options_.TeX && ~isoctave && matlab_ver_less_than('8.4') && ~isempty(getenv('TOP_TEST_DIR'))
    warning('Disabling TeX option due to a bug in MATLAB R2014a with -nodisplay')
    options_.TeX = false;
end

if M_.exo_nbr==0
    error('stoch_simul:: does not support having no varexo in the model. As a workaround you could define a dummy exogenous variable.')
end

if ismember(options_.solve_algo,[10,11])
    error('stoch_simul:: perturbation solutions are not compatible with mixed complementarity problems and their solvers')
end

test_for_deep_parameters_calibration(M_);

dr = oo_.dr;

options_old = options_;
if options_.linear
    options_.order = 1;
end
if options_.order == 1
    options_.replic = 1;
end

if options_.order~=1 && M_.hessian_eq_zero
    options_.order = 1;
    warning('stoch_simul: using order = 1 because Hessian is equal to zero');
end

if isempty(options_.qz_criterium)
    options_.qz_criterium = 1+1e-6;
end

if options_.partial_information || options_.ACES_solver
    PI_PCL_solver = 1;
    if options_.order ~= 1
        warning('stoch_simul:: forcing order=1 since you are using partial_information or ACES solver')
        options_.order = 1;
    end
else
    PI_PCL_solver = 0;
end

TeX = options_.TeX;

if isempty(var_list)
    var_list = M_.endo_names(1:M_.orig_endo_nbr);
end

[i_var, nvar, index_uniques] = varlist_indices(var_list, M_.endo_names);
var_list=var_list(index_uniques);
oo_.var_list = var_list;

iter_ = max(options_.periods,1);
if M_.exo_nbr > 0
    oo_.exo_simul= ones(iter_ + M_.maximum_lag + M_.maximum_lead,1) * oo_.exo_steady_state';
end

check_model(M_);

oo_.dr=set_state_space(dr,M_,options_);

if PI_PCL_solver
    [oo_.dr, info] = PCL_resol(oo_.steady_state,0);
elseif options_.discretionary_policy
    if ~options_.order==1
        error('discretionary_policy: only linear-quadratic problems can be solved');
    end
    [~,info,M_,oo_] = discretionary_policy_1(options_.instruments,M_,options_,oo_);
else
    if options_.logged_steady_state %if steady state was previously logged, undo this
        oo_.dr.ys=exp(oo_.dr.ys);
        oo_.steady_state=exp(oo_.steady_state);
        options_.logged_steady_state=0;
    end
    [~,info,M_,oo_] = resol(0,M_,options_,oo_);
end

if options_.loglinear && isfield(oo_.dr,'ys') && options_.logged_steady_state==0 %log steady state for correct display of decision rule
    oo_.dr.ys=log_variable(1:M_.endo_nbr,oo_.dr.ys,M_);
    oo_.steady_state=log_variable(1:M_.endo_nbr,oo_.steady_state,M_);
    options_old.logged_steady_state = 1; %make sure option is preserved outside of stoch_simul
    options_.logged_steady_state=1; %set option for use in stoch_simul
end
if info(1)
    options_ = options_old;
    print_info(info, options_.noprint, options_);
    return
end

if ~options_.noprint
    skipline()
    disp('MODEL SUMMARY')
    skipline()
    disp(['  Number of variables:         ' int2str(M_.endo_nbr)])
    disp(['  Number of stochastic shocks: ' int2str(M_.exo_nbr)])
    disp(['  Number of state variables:   ' int2str(M_.nspred)])
    disp(['  Number of jumpers:           ' int2str(M_.nsfwrd)])
    disp(['  Number of static variables:  ' int2str(M_.nstatic)])
    my_title='MATRIX OF COVARIANCE OF EXOGENOUS SHOCKS';
    labels = M_.exo_names;
    headers = vertcat('Variables', labels);
    lh = cellofchararraymaxlength(labels)+2;
    dyntable(options_, my_title, headers, labels, M_.Sigma_e, lh, 10, 6);
    if options_.TeX
        labels = M_.exo_names_tex;
        headers = vertcat('Variables', labels);
        lh = cellofchararraymaxlength(labels)+2;
        dyn_latex_table(M_, options_, my_title, 'covar_ex_shocks', headers, labels, M_.Sigma_e, lh, 10, 6);
    end
    if ~all(diag(M_.H)==0)
        my_title='MATRIX OF COVARIANCE OF MEASUREMENT ERRORS';
        labels = cellfun(@(x) horzcat('SE_', x), options_.varobs, 'UniformOutput', false);
        headers = vertcat('Variables', labels);
        lh = cellofchararraymaxlength(labels)+2;
        dyntable(options_, my_title, headers, labels, M_.H, lh, 10, 6);
        if options_.TeX
            labels = M_.exo_names_tex;
            headers = vertcat('Variables', labels);
            lh = cellofchararraymaxlength(labels)+2;
            dyn_latex_table(M_, options_, my_title, 'covar_ME', headers, labels, M_.H, lh, 10, 6);
        end
    end
    if options_.partial_information
        skipline()
        disp('SOLUTION UNDER PARTIAL INFORMATION')
        skipline()

        if isfield(options_,'varobs')&& ~isempty(options_.varobs)
            PCL_varobs = options_.varobs;
            disp('OBSERVED VARIABLES')
        else
            PCL_varobs = M_.endo_names;
            disp(' VAROBS LIST NOT SPECIFIED')
            disp(' ASSUMED OBSERVED VARIABLES')
        end
        for i=1:length(PCL_varobs)
            disp(['    ' PCL_varobs{i}])
        end
    end
    skipline()
    if options_.order <= 2 && ~PI_PCL_solver
        if ~options_.nofunctions
            disp_dr(oo_.dr,options_.order,var_list);
        end
    end
end

if options_.periods > 0 && ~PI_PCL_solver
    if options_.periods <= options_.drop
        fprintf('\nSTOCH_SIMUL error: The horizon of simulation is shorter than the number of observations to be dropped.\n')
        fprintf('STOCH_SIMUL error: Either increase options_.periods or decrease options_.drop.\n')
        options_ =options_old;
        return
    end
    if isempty(M_.endo_histval)
        y0 = oo_.dr.ys;
    else
        if options_.loglinear
            y0 = log_variable(1:M_.endo_nbr,M_.endo_histval,M_);
        else
            y0 = M_.endo_histval;
        end
    end
    [ys, oo_] = simult(y0,oo_.dr,M_,options_,oo_);
    oo_.endo_simul = ys;
    if ~options_.minimal_workspace
        dyn2vec(M_, oo_, options_);
    end
end

if ~options_.nomoments
    if PI_PCL_solver
        PCL_Part_info_moments(0, PCL_varobs, oo_.dr, i_var);
    elseif options_.periods == 0
        if options_.order == 1 || (options_.order == 2 && ~options_.pruning)
            oo_=disp_th_moments(oo_.dr,var_list,M_,options_,oo_);
        elseif (ismember(options_.order,[2,3])) && options_.pruning  
            % There is no code for theoretical moments at 3rd order without pruning
            oo_=disp_th_moments_pruned_state_space(oo_.dr,M_,options_,i_var,oo_);
        end
    else
        oo_=disp_moments(oo_.endo_simul,var_list,M_,options_,oo_);
    end
end


if options_.irf
    var_listTeX = M_.endo_names_tex(i_var);
    if ~options_.nograph || (TeX && any(strcmp('eps',cellstr(options_.graph_format))))
        if ~exist([M_.dname '/graphs'],'dir')
            mkdir(M_.dname,'graphs');
        end
    end
    if TeX && any(strcmp('eps',cellstr(options_.graph_format)))
        fidTeX = fopen([M_.dname, '/graphs/' M_.fname '_IRF.tex'],'w');
        fprintf(fidTeX,'%% TeX eps-loader file generated by stoch_simul.m (Dynare).\n');
        fprintf(fidTeX,['%% ' datestr(now,0) '\n']);
        fprintf(fidTeX,' \n');
    end
    SS(M_.exo_names_orig_ord,M_.exo_names_orig_ord)=M_.Sigma_e+1e-14*eye(M_.exo_nbr);
    cs = transpose(chol(SS));
    tit(M_.exo_names_orig_ord) = M_.exo_names;
    if TeX
        titTeX(M_.exo_names_orig_ord) = M_.exo_names_tex;
    end
    irf_shocks_indx = getIrfShocksIndx();
    for i=irf_shocks_indx
        if SS(i,i) > 1e-13
            if PI_PCL_solver
                y=PCL_Part_info_irf (0, PCL_varobs, i_var, M_, oo_.dr, options_.irf, i);
            else
                if options_.order>1 && options_.relative_irf % normalize shock to 0.01 before IRF generation for GIRFs; multiply with 100 later
                    y=irf(M_, options_, oo_.dr,cs(M_.exo_names_orig_ord,i)./cs(i,i)/100, options_.irf, options_.drop, ...
                          options_.replic, options_.order);
                else %for linear model, rescaling is done later
                    y=irf(M_, options_, oo_.dr,cs(M_.exo_names_orig_ord,i), options_.irf, options_.drop, ...
                          options_.replic, options_.order);
                end
            end
            if ~options_.noprint && any(any(isnan(y))) && ~options_.pruning && ~(options_.order==1)
                fprintf('\nstoch_simul:: The simulations conducted for generating IRFs to %s were explosive.\n', M_.exo_names{i})
                fprintf('stoch_simul:: No IRFs will be displayed. Either reduce the shock size, \n')
                fprintf('stoch_simul:: use pruning, or set the approximation order to 1.');
                skipline(2);
            end
            if options_.relative_irf
                if options_.order==1 %multiply with 100 for backward compatibility
                    y = 100*y/cs(i,i);
                end
            end
            irfs   = [];
            mylist = [];
            if TeX
                mylistTeX = [];
            end
            for j = 1:nvar
                assignin('base',[M_.endo_names{i_var(j)} '_' M_.exo_names{i}],...
                         y(i_var(j),:)');
                oo_.irfs.([M_.endo_names{i_var(j)} '_' M_.exo_names{i}]) = y(i_var(j),:);
                if max(abs(y(i_var(j),:))) >= options_.impulse_responses.plot_threshold
                    irfs  = cat(1,irfs,y(i_var(j),:));
                    if isempty(mylist)
                        mylist = var_list{j};
                    else
                        mylist = char(mylist, var_list{j});
                    end
                    if TeX
                        if isempty(mylistTeX)
                            mylistTeX = var_listTeX{j};
                        else
                            mylistTeX = char(mylistTeX, var_listTeX{j});
                        end
                    end
                else
                    if options_.debug
                        fprintf('stoch_simul:: The IRF of %s to %s is smaller than the irf_plot_threshold of %4.3f and will not be displayed.\n',M_.endo_names{i_var(j)},M_.exo_names{i},options_.impulse_responses.plot_threshold)
                    end
                end
            end
            if ~options_.nograph
                number_of_plots_to_draw = size(irfs,1);
                [nbplt,nr,nc,lr,lc,nstar] = pltorg(number_of_plots_to_draw);
                if nbplt == 0
                elseif nbplt == 1
                    if options_.relative_irf
                        hh = dyn_figure(options_.nodisplay,'Name',['Relative response to' ...
                                            ' orthogonalized shock to ' tit{i}]);
                    else
                        hh = dyn_figure(options_.nodisplay,'Name',['Orthogonalized shock to' ...
                                            ' ' tit{i}]);
                    end
                    for j = 1:number_of_plots_to_draw
                        subplot(nr,nc,j);
                        plot(1:options_.irf,transpose(irfs(j,:)),'-k','linewidth',1);
                        hold on
                        plot([1 options_.irf],[0 0],'-r','linewidth',0.5);
                        hold off
                        xlim([1 options_.irf]);
                        remove_fractional_xticks;
                        if TeX
                            title(['$' deblank(mylistTeX(j,:)) '$'],'Interpreter','latex');
                        else
                            title(deblank(mylist(j,:)),'Interpreter','none');
                        end
                    end
                    dyn_saveas(hh,[M_.dname, '/graphs/' M_.fname '_IRF_' tit{i}],options_.nodisplay,options_.graph_format);
                    if TeX && any(strcmp('eps',cellstr(options_.graph_format)))
                        fprintf(fidTeX,'\\begin{figure}[H]\n');
                        fprintf(fidTeX,'\\centering \n');
                        fprintf(fidTeX,'\\includegraphics[width=%2.2f\\textwidth]{%s_IRF_%s}\n',options_.figures.textwidth*min(j/nc,1),[M_.dname, '/graphs/' M_.fname],tit{i});
                        fprintf(fidTeX,'\\caption{Impulse response functions (orthogonalized shock to $%s$).}\n',titTeX{i});
                        fprintf(fidTeX,'\\label{Fig:IRF:%s}\n', tit{i});
                        fprintf(fidTeX,'\\end{figure}\n');
                        fprintf(fidTeX,' \n');
                    end
                else
                    for fig = 1:nbplt-1
                        if options_.relative_irf
                            hh = dyn_figure(options_.nodisplay,'Name',['Relative response to orthogonalized shock' ...
                                                ' to ' tit{i} ' figure ' int2str(fig)]);
                        else
                            hh = dyn_figure(options_.nodisplay,'Name',['Orthogonalized shock to ' tit{i} ...
                                                ' figure ' int2str(fig)]);
                        end
                        for plt = 1:nstar
                            subplot(nr,nc,plt);
                            plot(1:options_.irf,transpose(irfs((fig-1)*nstar+plt,:)),'-k','linewidth',1);
                            hold on
                            plot([1 options_.irf],[0 0],'-r','linewidth',0.5);
                            hold off
                            xlim([1 options_.irf]);
                            remove_fractional_xticks
                            if TeX
                                title(['$' deblank(mylistTeX((fig-1)*nstar+plt,:)) '$'],'Interpreter','latex');
                            else
                                title(deblank(mylist((fig-1)*nstar+plt,:)),'Interpreter','none');
                            end
                        end
                        dyn_saveas(hh,[M_.dname, '/graphs/'  M_.fname '_IRF_' tit{i} int2str(fig)],options_.nodisplay,options_.graph_format);
                        if TeX && any(strcmp('eps',cellstr(options_.graph_format)))
                            fprintf(fidTeX,'\\begin{figure}[H]\n');
                            fprintf(fidTeX,'\\centering \n');
                            fprintf(fidTeX,'\\includegraphics[width=%2.2f\\textwidth]{%s_IRF_%s%s}\n',options_.figures.textwidth*min(plt/nc,1),[M_.dname, '/graphs/' M_.fname],tit{i},int2str(fig));
                            if options_.relative_irf
                                fprintf(fidTeX,'\\caption{Relative impulse response functions (orthogonalized shock to $%s$).}', titTeX{i});
                            else
                                fprintf(fidTeX,'\\caption{Impulse response functions (orthogonalized shock to $%s$).}', titTeX{i});
                            end
                            fprintf(fidTeX,'\\label{Fig:IRF:%s:%s}\n', tit{i},int2str(fig));
                            fprintf(fidTeX,'\\end{figure}\n');
                            fprintf(fidTeX,' \n');
                        end
                    end
                    hh = dyn_figure(options_.nodisplay,'Name',['Orthogonalized shock to ' tit{i} ' figure ' int2str(nbplt) '.']);
                    m = 0;
                    for plt = 1:number_of_plots_to_draw-(nbplt-1)*nstar
                        m = m+1;
                        subplot(lr,lc,m);
                        plot(1:options_.irf,transpose(irfs((nbplt-1)*nstar+plt,:)),'-k','linewidth',1);
                        hold on
                        plot([1 options_.irf],[0 0],'-r','linewidth',0.5);
                        hold off
                        xlim([1 options_.irf]);
                        remove_fractional_xticks
                            if TeX
                                title(['$' deblank(mylistTeX((nbplt-1)*nstar+plt,:)) '$'],'Interpreter','latex');
                            else
                                title(deblank(mylist((nbplt-1)*nstar+plt,:)),'Interpreter','none');
                            end
                    end
                    dyn_saveas(hh,[M_.dname, '/graphs/' M_.fname '_IRF_' tit{i} int2str(nbplt) ],options_.nodisplay,options_.graph_format);
                    if TeX && any(strcmp('eps',cellstr(options_.graph_format)))
                        fprintf(fidTeX,'\\begin{figure}[H]\n');
                        fprintf(fidTeX,'\\centering \n');
                        fprintf(fidTeX,'\\includegraphics[width=%2.2f\\textwidth]{%s_IRF_%s%s}\n',options_.figures.textwidth*min(m/lc,1),[M_.dname, '/graphs/' M_.fname],tit{i},int2str(nbplt));
                        if options_.relative_irf
                            fprintf(fidTeX,'\\caption{Relative impulse response functions (orthogonalized shock to $%s$).}', titTeX{i});
                        else
                            fprintf(fidTeX,'\\caption{Impulse response functions (orthogonalized shock to $%s$).}', titTeX{i});
                        end
                        fprintf(fidTeX,'\\label{Fig:IRF:%s:%s}\n', tit{i},int2str(nbplt));
                        fprintf(fidTeX,'\\end{figure}\n');
                        fprintf(fidTeX,' \n');
                    end
                end
            end
        end
    end
    if TeX && any(strcmp('eps',cellstr(options_.graph_format)))
        fprintf(fidTeX,' \n');
        fprintf(fidTeX,'%% End Of TeX file. \n');
        fclose(fidTeX);
    end
end

if options_.SpectralDensity.trigger
    [oo_] = UnivariateSpectralDensity(M_,oo_,options_,var_list);
end


options_ = options_old;
% temporary fix waiting for local options
options_.partial_information = 0;

oo_.gui.ran_stoch_simul = true;
