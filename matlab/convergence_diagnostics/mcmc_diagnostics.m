function oo_ = mcmc_diagnostics(options_, estim_params_, M_, oo_)
% function oo_ = mcmc_diagnostics(options_, estim_params_, M_, oo_)
% Computes convergence tests
%
% INPUTS
%   options_         [structure]    Dynare options structure
%   estim_params_    [structure]    Dynare estimation parameter structure
%   M_               [structure]    Dynare model structure
%   oo_              [structure]    Dynare results structure
%
% OUTPUTS
%   oo_              [structure]
%
% SPECIAL REQUIREMENTS
%   none
%
% PARALLEL CONTEXT
% See the comment in posterior_sampler.m funtion.

% Copyright © 2005-2023 Dynare Team
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
graphFolder = CheckPath('graphs',M_.dname);
latexFolder = CheckPath('latex',M_.dname);
MetropolisFolder = CheckPath('metropolis',M_.dname);
ModelName = M_.fname;

TeX = options_.TeX;

record=load_last_mh_history_file(MetropolisFolder, ModelName);
NumberOfMcFilesPerBlock = record.LastFileNumber;
[nblck, npar] = size(record.LastParameters);
npardisp = options_.convergence.brooksgelman.plotrows;

% check if all the mh files are available.
issue_an_error_message = 0;
for b = 1:nblck
    nfiles = length(dir([MetropolisFolder ,filesep, ModelName '_mh*_blck' num2str(b) '.mat']));
    if ~isequal(NumberOfMcFilesPerBlock,nfiles)
        issue_an_error_message = 1;
        fprintf('The number of MCMC files in chain %u is %u while the mh-history files indicate that we should have %u MCMC files per chain!\n',b, nfiles, NumberOfMcFilesPerBlock);
    end
end
if issue_an_error_message
    error('mcmc_diagnostics: I cannot proceed because some MCMC files are missing. Check your MCMC files...!');
end

% compute inefficiency factor
FirstLine = record.KeepedDraws.FirstLine;
TotalNumberOfMhFiles = sum(record.MhDraws(:,2));
TotalNumberOfMhDraws = sum(record.MhDraws(:,1));
FirstMhFile = record.KeepedDraws.FirstMhFile;
NumberOfDraws = TotalNumberOfMhDraws-floor(options_.mh_drop*TotalNumberOfMhDraws);

param_name = {};
param_name_tex = {};

Ifac=NaN(nblck,npar);
for jj = 1:npar
    if options_.TeX
        [par_name_temp, par_name_tex_temp] = get_the_name(jj, options_.TeX, M_,estim_params_, options_.varobs);
        param_name = vertcat(param_name, par_name_temp);
        par_name_tex_temp = strrep(par_name_tex_temp,'$','');
        param_name_tex = vertcat(param_name_tex, par_name_tex_temp);
    else
        par_name_temp = get_the_name(jj, options_.TeX, M_, estim_params_, options_.varobs);
        param_name = vertcat(param_name, par_name_temp);
    end
    Draws = GetAllPosteriorDraws(options_, M_.dname, M_.fname, jj, FirstMhFile, FirstLine, TotalNumberOfMhFiles, NumberOfDraws, nblck);
    Draws = reshape(Draws, [NumberOfDraws nblck]);
    Nc = min(1000, NumberOfDraws/2);
    for ll = 1:nblck
        Ifac(ll,jj) = mcmc_ifac(Draws(:,ll), Nc);
    end
    tmp = num2cell(Ifac(:,jj));
end

my_title='MCMC Inefficiency factors per block';
IFAC_header = {'Parameter'};
IFAC_header_tex = {'Parameter'};
for j = 1:nblck
    IFAC_header = vertcat(IFAC_header, ['Block ' int2str(j)]);
    IFAC_header_tex = vertcat(IFAC_header_tex, ['Block~' int2str(j)]);
end

lh = cellofchararraymaxlength(param_name)+2;
dyntable(options_, my_title, IFAC_header, param_name, Ifac', lh, 12, 3);
skipline()

if options_.TeX
    dyn_latex_table(M_, options_, my_title, 'MCMC_inefficiency_factors', IFAC_header_tex, param_name_tex, Ifac', lh, 12, 3);
end
record.InefficiencyFactorsPerBlock = Ifac;
update_last_mh_history_file(MetropolisFolder, ModelName, record);


PastDraws = sum(record.MhDraws,1);
NumberOfDraws  = PastDraws(1);

if ~strcmp(options_.posterior_sampler_options.posterior_sampling_method,'slice') && NumberOfDraws<=2000
    warning('MCMC convergence diagnostics are not computed because the total number of iterations is not bigger than 2000!');
    return
end

convergence_diagnostics_geweke=zeros(npar,4+2*length(options_.convergence.geweke.taper_steps));
if any(options_.convergence.geweke.geweke_interval<0) || any(options_.convergence.geweke.geweke_interval>1) || length(options_.convergence.geweke.geweke_interval)~=2 ...
        || (options_.convergence.geweke.geweke_interval(2)-options_.convergence.geweke.geweke_interval(1)<0)
    fprintf('\nCONVERGENCE DIAGNOSTICS: Invalid option for geweke_interval. Using the default of [0.2 0.5].\n')
    options_.convergence.geweke.geweke_interval=[0.2 0.5];
end
first_obs_begin_sample = max(1,ceil(options_.mh_drop*NumberOfDraws));
last_obs_begin_sample = first_obs_begin_sample+round(options_.convergence.geweke.geweke_interval(1)*NumberOfDraws*(1-options_.mh_drop));
first_obs_end_sample = first_obs_begin_sample+round(options_.convergence.geweke.geweke_interval(2)*NumberOfDraws*(1-options_.mh_drop));
param_name = {};
if options_.TeX
    param_name_tex = {};
end
for jj=1:npar
    if options_.TeX
        [param_name_temp, param_name_tex_temp] = get_the_name(jj, options_.TeX, M_, estim_params_, options_.varobs);
        param_name_tex = vertcat(param_name_tex, strrep(param_name_tex_temp, '$',''));
        param_name = vertcat(param_name, param_name_temp);
    else
        param_name_temp = get_the_name(jj, options_.TeX, M_,estim_params_, options_.varobs);
        param_name = vertcat(param_name, param_name_temp);
    end
end
datamat=NaN(npar,3+length(options_.convergence.geweke.taper_steps),nblck);
%remove stale results as it will cause assignment problems
if options_.convergence.rafterylewis.indicator && isfield(oo_,'convergence') && isfield(oo_.convergence,'raftery_lewis')
    oo_.convergence=rmfield(oo_.convergence,'raftery_lewis');
end

for block_iter=1:nblck
    fprintf('\n\nConvergence diagnostics results for chain %u.\n',block_iter);
    fprintf('\nGeweke (1992) Convergence Tests, based on means of draws %d to %d vs %d to %d for chain %u.\n',first_obs_begin_sample,last_obs_begin_sample,first_obs_end_sample,NumberOfDraws,block_iter);
    fprintf('p-values are for Chi2-test for equality of means.\n');
    Geweke_header = {'Parameter'; 'Post. Mean'; 'Post. Std'; 'p-val No Taper'};
    for ii = 1:length(options_.convergence.geweke.taper_steps)
        Geweke_header = vertcat(Geweke_header, ['p-val ' num2str(options_.convergence.geweke.taper_steps(ii)),'% Taper']);
    end
    for jj=1:npar
        startline=0;
        for n = 1:NumberOfMcFilesPerBlock
            load([MetropolisFolder '/' ModelName '_mh',int2str(n),'_blck',num2str(block_iter),'.mat'],'x2');
            nx2 = size(x2,1);
            param_draws(startline+(1:nx2),1) = x2(:,jj);
            startline = startline + nx2;
        end
        [results_vec, results_struct] = geweke_moments(param_draws,options_);
        convergence_diagnostics_geweke(jj,:)=results_vec;

        param_draws1 = param_draws(first_obs_begin_sample:last_obs_begin_sample,:);
        param_draws2 = param_draws(first_obs_end_sample:end,:);
        [results_vec1] = geweke_moments(param_draws1,options_);
        [results_vec2] = geweke_moments(param_draws2,options_);

        results_struct = geweke_chi2_test(results_vec1,results_vec2,results_struct,options_);
        oo_.convergence.geweke(block_iter).(param_name{jj}) = results_struct;
        datamat(jj,:,block_iter)=[results_struct.posteriormean,results_struct.posteriorstd,results_struct.prob_chi2_test];
    end
    lh = size(param_name,2)+2;
    dyntable(options_, '', Geweke_header, param_name, datamat(:,:,block_iter), lh, 12, 3);
    if options_.TeX
        Geweke_tex_header = {'Parameter'; 'Mean'; 'Std'; 'No\ Taper'};
        additional_header = {[' & \multicolumn{2}{c}{Posterior} & \multicolumn{',num2str(1+length(options_.convergence.geweke.taper_steps)),'}{c}{p-values} \\'],
            ['\cmidrule(r{.75em}){2-3} \cmidrule(r{.75em}){4-',num2str(4+length(options_.convergence.geweke.taper_steps)),'}']};
        for ii=1:length(options_.convergence.geweke.taper_steps)
            Geweke_tex_header = vertcat(Geweke_tex_header, [num2str(options_.convergence.geweke.taper_steps(ii)),'\%%\ Taper']);
        end
        headers = Geweke_tex_header;
        lh = cellofchararraymaxlength(param_name_tex)+2;
        my_title=sprintf('Geweke (1992) Convergence Tests, based on means of draws %d to %d vs %d to %d for chain %u. p-values are for $\\\\chi^2$-test for equality of means.',first_obs_begin_sample,last_obs_begin_sample,first_obs_end_sample,NumberOfDraws,block_iter);
        dyn_latex_table(M_, options_, my_title, ['geweke_block_' num2str(block_iter)], headers, param_name_tex, datamat(:,:,block_iter), lh, 12, 4, additional_header);
    end
    skipline(2);

    if options_.convergence.rafterylewis.indicator
        if any(options_.convergence.rafterylewis.qrs<0) || any(options_.convergence.rafterylewis.qrs>1) || length(options_.convergence.rafterylewis.qrs)~=3 ...
                || (options_.convergence.rafterylewis.qrs(1)-options_.convergence.rafterylewis.qrs(2)<=0)
            fprintf('\nInvalid option for raftery_lewis_qrs. Using the default of [0.025 0.005 0.95].\n');
            options_.convergence.rafterylewis.qrs=[0.025 0.005 0.95];
        end
        Raftery_Lewis_q=options_.convergence.rafterylewis.qrs(1);
        Raftery_Lewis_r=options_.convergence.rafterylewis.qrs(2);
        Raftery_Lewis_s=options_.convergence.rafterylewis.qrs(3);
        oo_.convergence.raftery_lewis(block_iter) = raftery_lewis(x2,Raftery_Lewis_q,Raftery_Lewis_r,Raftery_Lewis_s);
        my_title=sprintf('Raftery/Lewis (1992) Convergence Diagnostics, based on quantile q=%4.3f with precision r=%4.3f with probability s=%4.3f for chain %u.',Raftery_Lewis_q,Raftery_Lewis_r,Raftery_Lewis_s,block_iter);
        headers = {'Variables'; 'M (burn-in)'; 'N (req. draws)'; 'N+M (total draws)'; 'k (thinning)'};
        raftery_data_mat=[oo_.convergence.raftery_lewis(block_iter).M_burn,oo_.convergence.raftery_lewis(block_iter).N_prec,oo_.convergence.raftery_lewis(block_iter).N_total,oo_.convergence.raftery_lewis(block_iter).k_thin];
        raftery_data_mat=[raftery_data_mat; max(raftery_data_mat,[],1)];
        labels_Raftery_Lewis = vertcat(param_name, 'Maximum');
        lh = cellofchararraymaxlength(labels_Raftery_Lewis)+2;
        dyntable(options_, my_title, headers, labels_Raftery_Lewis, raftery_data_mat, lh, 10, 0);
        if options_.TeX
            labels_Raftery_Lewis_tex = vertcat(param_name_tex, 'Maximum');
            lh = cellofchararraymaxlength(labels_Raftery_Lewis_tex)+2;
            dyn_latex_table(M_, options_, my_title, ['raftery_lewis_' num2str(block_iter)], headers, labels_Raftery_Lewis_tex, raftery_data_mat, lh, 10, 0);
        end
    end
end
for block_iter=1:nblck
    oo_.convergence.raftery_lewis(block_iter).parameter_names=param_name;
end
if nblck==1
    return
end

if strcmp(options_.posterior_sampler_options.posterior_sampling_method,'slice') && NumberOfDraws<2000
    Origin = 1;
    StepSize = 1;
else
    Origin = 1000;
    StepSize = ceil((NumberOfDraws-Origin)/100);% So that the computational time does not
end
ALPHA = 0.2;                                % increase too much with the number of simulations.
time = 1:NumberOfDraws;
xx = Origin:StepSize:NumberOfDraws;
NumberOfLines = length(xx);

if NumberOfDraws < Origin
    fprintf('The number of simulations is too small to compute the MCMC convergence diagnostics.\n');
    return
end

if TeX && any(strcmp('eps',cellstr(options_.graph_format)))
    fidTeX = fopen([latexFolder '/' ModelName '_UnivariateDiagnostics.tex'],'w');
    fprintf(fidTeX,'%% TeX eps-loader file generated by mcmc_diagnostics.m (Dynare).\n');
    fprintf(fidTeX,['%% ' datestr(now,0) '\n']);
    fprintf(fidTeX,' \n');
end

fprintf('Univariate convergence diagnostic, Brooks and Gelman (1998):\n');

% The mandatory variables for local/remote parallel
% computing are stored in localVars struct.

localVars.MetropolisFolder = MetropolisFolder;
localVars.nblck = nblck;
localVars.NumberOfMcFilesPerBlock = NumberOfMcFilesPerBlock;
localVars.Origin = Origin;
localVars.StepSize = StepSize;
localVars.mh_drop = options_.mh_drop;
localVars.NumberOfDraws = NumberOfDraws;
localVars.NumberOfLines = NumberOfLines;
localVars.time = time;
localVars.M_ = M_;


% Like sequential execution!
if isnumeric(options_.parallel)
    fout = mcmc_diagnostics_core(localVars,1,npar,0);
    UDIAG = fout.UDIAG;
    clear fout
    % Parallel execution!
else
    if ~isempty(M_.bvar)
        ModelName = [ModelName '_bvar'];
    end
    NamFileInput={[M_.dname '/metropolis/'],[ModelName '_mh*_blck*.mat']};

    [fout, nBlockPerCPU, totCPU] = masterParallel(options_.parallel, 1, npar,NamFileInput,'mcmc_diagnostics_core', localVars, [], options_.parallel_info);
    UDIAG = fout(1).UDIAG;
    for j=2:totCPU
        UDIAG = cat(3,UDIAG ,fout(j).UDIAG);
    end
end

UDIAG(:,[2 4 6],:) = UDIAG(:,[2 4 6],:)/nblck;
skipline()
clear pmet temp moyenne CSUP CINF csup cinf n linea iter tmp;

pages = floor(npar/npardisp); % changed from 3 to npardisp
k = 0;
for i = 1:pages
    hh_fig = dyn_figure(options_.nodisplay,'Name','MCMC univariate convergence diagnostic (Brooks and Gelman,1998)');
    boxplot = 1;
    for j = 1:npardisp % Loop over parameters  %npardisp instead of 3
        k = k+1;
        [nam,namtex] = get_the_name(k,TeX,M_,estim_params_,options_.varobs);
        for crit = 1:3% Loop over criteria
            if crit == 1
                plt1 = UDIAG(:,1,k);
                plt2 = UDIAG(:,2,k);
                namnam  = [nam , ' (Interval)'];
            elseif crit == 2
                plt1 = UDIAG(:,3,k);
                plt2 = UDIAG(:,4,k);
                namnam  = [nam , ' (m2)'];
            elseif crit == 3
                plt1 = UDIAG(:,5,k);
                plt2 = UDIAG(:,6,k);
                namnam  = [nam , ' (m3)'];
            end
            subplot(npardisp,3,boxplot);  %Added more rows to display more variables
            plot(xx,plt1,'-b');     % Pooled
            hold on;
            plot(xx,plt2,'-r');     % Within (mean)
            hold off;
            xlim([xx(1) xx(NumberOfLines)])
            if TeX
                title(namtex,'interpreter','latex')
            else
                title(namnam,'Interpreter','none')
            end

            boxplot = boxplot + 1;
        end
    end
    dyn_saveas(hh_fig,[graphFolder '/' ModelName '_udiag' int2str(i)],options_.nodisplay,options_.graph_format);
    if TeX && any(strcmp('eps',cellstr(options_.graph_format)))
        fprintf(fidTeX,'\\begin{figure}[H]\n');
        fprintf(fidTeX,'\\centering \n');
        fprintf(fidTeX,'\\includegraphics[width=%2.2f\\textwidth]{%s_udiag%s}\n',options_.figures.textwidth*min((boxplot-1)/3,1),[graphFolder '/' ModelName],int2str(i));
        fprintf(fidTeX,'\\caption{Univariate convergence diagnostics for the Metropolis-Hastings.\n');
        fprintf(fidTeX,'The first, second and third columns are respectively the criteria based on\n');
        fprintf(fidTeX,'the eighty percent interval, the second and third moments.}');
        fprintf(fidTeX,'\\label{Fig:UnivariateDiagnostics:%s}\n',int2str(i));
        fprintf(fidTeX,'\\end{figure}\n');
        fprintf(fidTeX,'\n');
    end
end
reste = npar-k;
if reste
    if reste == 1
        nr = 3;
        nc = 1;
    else % Conditional for additional rows (variables) when not towards the end of the loop
        nr = npardisp;
        nc = 3;
    end
    hh_fig = dyn_figure(options_.nodisplay,'Name','MCMC univariate convergence diagnostic (Brooks and Gelman, 1998)');
    boxplot = 1;
    for j = 1:reste
        k = k+1;
        [nam,namtex] = get_the_name(k,TeX,M_,estim_params_,options_.varobs);
        for crit = 1:3
            if crit == 1
                plt1 = UDIAG(:,1,k);
                plt2 = UDIAG(:,2,k);
                namnam  = [nam , ' (Interval)'];
                if TeX
                    namnamtex  = [namtex , ' (Interval)'];
                end
            elseif crit == 2
                plt1 = UDIAG(:,3,k);
                plt2 = UDIAG(:,4,k);
                namnam  = [nam , ' (m2)'];
                if TeX
                    namnamtex  = [namtex , ' (m2)'];
                end
            elseif crit == 3
                plt1 = UDIAG(:,5,k);
                plt2 = UDIAG(:,6,k);
                namnam  = [nam , ' (m3)'];
                if TeX
                    namnamtex  = [namtex , ' (m3)'];
                end
            end
            subplot(nr,nc,boxplot);
            plot(xx,plt1,'-b');                                 % Pooled
            hold on;
            plot(xx,plt2,'-r');                                 % Within (mean)
            hold off;
            xlim([xx(1) xx(NumberOfLines)]);
            if TeX
                title(namnamtex,'Interpreter','latex');
            else
                title(namnam,'Interpreter','none');
            end
            boxplot = boxplot + 1;
        end
    end
    dyn_saveas(hh_fig,[ graphFolder '/' ModelName '_udiag' int2str(pages+1)],options_.nodisplay,options_.graph_format);
    if TeX && any(strcmp('eps',cellstr(options_.graph_format)))
        fprintf(fidTeX,'\\begin{figure}[H]\n');
        fprintf(fidTeX,'\\centering \n');
        fprintf(fidTeX,'\\includegraphics[width=%2.2f\\textwidth]{%s_udiag%s}\n',options_.figures.textwidth*min((boxplot-1)/nc,1),[graphFolder '/' ModelName],int2str(pages+1));
        if reste == 2
            fprintf(fidTeX,'\\caption{Univariate convergence diagnostics for the Metropolis-Hastings.\n');
            fprintf(fidTeX,'The first, second and third columns are respectively the criteria based on\n');
            fprintf(fidTeX,'the eighty percent interval, the second and third moments.}');
        elseif reste == 1
            fprintf(fidTeX,'\\caption{Univariate convergence diagnostics for the Metropolis-Hastings.\n');
            fprintf(fidTeX,'The first, second and third rows are respectively the criteria based on\n');
            fprintf(fidTeX,'the eighty percent interval, the second and third moments.}');
        end
        fprintf(fidTeX,'\\label{Fig:UnivariateDiagnostics:%s}\n',int2str(pages+1));
        fprintf(fidTeX,'\\end{figure}\n');
        fprintf(fidTeX,'\n');
        fprintf(fidTeX,'%% End Of TeX file.');
        fclose(fidTeX);
    end
end % if reste > 0
clear UDIAG;
%
% Multivariate diagnostic.
%
if TeX && any(strcmp('eps',cellstr(options_.graph_format)))
    fidTeX = fopen([latexFolder '/' ModelName '_MultivariateDiagnostics.tex'],'w');
    fprintf(fidTeX,'%% TeX eps-loader file generated by mcmc_diagnostics.m (Dynare).\n');
    fprintf(fidTeX,['%% ' datestr(now,0) '\n']);
    fprintf(fidTeX,' \n');
end
tmp = zeros(NumberOfDraws*nblck,3);
MDIAG = zeros(NumberOfLines,6);
for b = 1:nblck
    startline = 0;
    for n = 1:NumberOfMcFilesPerBlock
        load([MetropolisFolder '/' ModelName '_mh',int2str(n),'_blck' int2str(b) '.mat'],'logpo2');
        nlogpo2 = size(logpo2,1);
        tmp((b-1)*NumberOfDraws+startline+(1:nlogpo2),1) = logpo2;
        startline = startline+nlogpo2;
    end
end
clear logpo2;
tmp(:,2) = kron(transpose(1:nblck),ones(NumberOfDraws,1));
tmp(:,3) = kron(ones(nblck,1),time');
tmp = sortrows(tmp,1);
ligne   = 0;
for iter  = Origin:StepSize:NumberOfDraws
    ligne = ligne+1;
    linea = ceil(options_.mh_drop*iter);
    n     = iter-linea+1;
    cinf  = max(1,round(n*ALPHA/2));
    csup  = round(n*(1-ALPHA/2));
    CINF  = max(1,round(nblck*n*ALPHA/2));
    CSUP  = round(nblck*n*(1-ALPHA/2));
    temp  = tmp(find((tmp(:,3)>=linea) & (tmp(:,3)<=iter)),1:2);
    MDIAG(ligne,1) = temp(CSUP,1)-temp(CINF,1);
    moyenne = mean(temp(:,1));%% Pooled mean.
    MDIAG(ligne,3) = sum((temp(:,1)-moyenne).^2)/(nblck*n-1);
    MDIAG(ligne,5) = sum(abs(temp(:,1)-moyenne).^3)/(nblck*n-1);
    for i=1:nblck
        pmet = temp(find(temp(:,2)==i));
        MDIAG(ligne,2) = MDIAG(ligne,2) + pmet(csup,1)-pmet(cinf,1);
        moyenne = mean(pmet,1); %% Within mean.
        MDIAG(ligne,4) = MDIAG(ligne,4) + sum((pmet(:,1)-moyenne).^2)/(n-1);
        MDIAG(ligne,6) = MDIAG(ligne,6) + sum(abs(pmet(:,1)-moyenne).^3)/(n-1);
    end
end
MDIAG(:,[2 4 6],:) = MDIAG(:,[2 4 6],:)/nblck;

hh_fig = dyn_figure(options_.nodisplay,'Name','Multivariate convergence diagnostic');
boxplot = 1;
for crit = 1:3
    if crit == 1
        plt1 = MDIAG(:,1);
        plt2 = MDIAG(:,2);
        namnam  = 'Interval';
    elseif crit == 2
        plt1 = MDIAG(:,3);
        plt2 = MDIAG(:,4);
        namnam  = 'm2';
    elseif crit == 3
        plt1 = MDIAG(:,5);
        plt2 = MDIAG(:,6);
        namnam  = 'm3';
    end
    subplot(3,1,boxplot);
    plot(xx,plt1,'-b');  % Pooled
    hold on
    plot(xx,plt2,'-r');  % Within (mean)
    hold off
    xlim([xx(1) xx(NumberOfLines)])
    title(namnam,'Interpreter','none');
    boxplot = boxplot + 1;
end
dyn_saveas(hh_fig,[ graphFolder '/' ModelName '_mdiag'],options_.nodisplay,options_.graph_format);

if TeX && any(strcmp('eps',cellstr(options_.graph_format)))
    fprintf(fidTeX,'\\begin{figure}[H]\n');
    fprintf(fidTeX,'\\centering \n');
    fprintf(fidTeX,'\\includegraphics[width=0.8\\textwidth]{%s_mdiag}\n',[graphFolder '/' ModelName]);
    fprintf(fidTeX,'\\caption{Multivariate convergence diagnostics for the Metropolis-Hastings.\n');
    fprintf(fidTeX,'The first, second and third rows are respectively the criteria based on\n');
    fprintf(fidTeX,'the eighty percent interval, the second and third moments. The different \n');
    fprintf(fidTeX,'parameters are aggregated using the posterior kernel.}');
    fprintf(fidTeX,'\\label{Fig:MultivariateDiagnostics}\n');
    fprintf(fidTeX,'\\end{figure}\n');
    fprintf(fidTeX,'\n');
    fprintf(fidTeX,'%% End Of TeX file.');
    fclose(fidTeX);
end
