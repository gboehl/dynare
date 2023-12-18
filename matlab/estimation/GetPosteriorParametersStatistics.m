function oo_ = GetPosteriorParametersStatistics(estim_params_, M_, options_, bayestopt_, oo_, pnames)

% This function prints and saves posterior estimates after the mcmc
% (+updates of oo_ & TeX output).
%
% INPUTS
%   estim_params_    [structure]    Dynare estimation parameter structure
%   M_               [structure]    Dynare model structure
%   options_         [structure]    Dynare options structure
%   bayestopt_       [structure]    Dynare structure describing priors
%   oo_              [structure]    Dynare results structure
%   pnames           [cell]         cell of strings, names of the prior shapes available
%
% OUTPUTS
%   oo_              [structure]    Dynare results structure
%
% SPECIAL REQUIREMENTS
%   None.

% Copyright © 2006-2023 Dynare Team
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

TeX     = options_.TeX;
nvx     = estim_params_.nvx;
nvn     = estim_params_.nvn;
ncx     = estim_params_.ncx;
ncn     = estim_params_.ncn;
np      = estim_params_.np ;

latexFolder = CheckPath('latex',M_.dname);
FileName = M_.fname;

if ~issmc(options_)
    MetropolisFolder = CheckPath('metropolis',M_.dname);
    record=load_last_mh_history_file(MetropolisFolder,FileName);
    FirstLine = record.KeepedDraws.FirstLine;
    TotalNumberOfMhFiles = sum(record.MhDraws(:,2));
    TotalNumberOfMhDraws = sum(record.MhDraws(:,1));
    FirstMhFile = record.KeepedDraws.FirstMhFile;
    NumberOfDraws = TotalNumberOfMhDraws-floor(options_.mh_drop*TotalNumberOfMhDraws);
    mh_nblck = size(record.LastParameters,1);
    clear record;
end

header_width = row_header_width(M_, estim_params_, bayestopt_);
hpd_interval=[num2str(options_.mh_conf_sig*100), '% HPD interval'];
tit2 = sprintf('%-*s %12s %12s %23s %8s %12s\n',header_width,' ','prior mean','post. mean',hpd_interval,'prior','pstdev');
pformat =      '%-*s %12.3f % 12.4f %11.4f %11.4f %8s %12.4f';

skipline(2)
disp('ESTIMATION RESULTS')
skipline()

if ishssmc(options_)
    dprintf('Log data density is %f.', oo_.MarginalDensity.hssmc);
    % Set function handle for GetAllPosteriorDraws
    getalldraws = @(i) GetAllPosteriorDraws(options_, M_.dname, [], i);
else
    if ~isfield(oo_,'MarginalDensity') || (issmc(options_) && ~isfield(oo_.MarginalDensity,'ModifiedHarmonicMean'))
        [~, oo_] = marginal_density(M_, options_, estim_params_, oo_, bayestopt_);
    end
    fprintf('Log data density is %f.', oo_.MarginalDensity.ModifiedHarmonicMean);
    % Set function handle for GetAllPosteriordraws
    getalldraws = @(i) GetAllPosteriorDraws(options_, M_.dname, M_.fname, i, FirstMhFile, FirstLine, TotalNumberOfMhFiles, NumberOfDraws, mh_nblck);
end

if ishssmc(options_)
    num_draws = options_.posterior_sampler_options.hssmc.particles;
    hpd_draws = round((1-options_.mh_conf_sig)*num_draws);
else
    num_draws=NumberOfDraws*mh_nblck;
    hpd_draws = round((1-options_.mh_conf_sig)*num_draws);
end

if hpd_draws<2
    fprintf('posterior_moments: There are not enough draws computes to compute HPD Intervals. Skipping their computation.\n')
end
if num_draws<9
    fprintf('posterior_moments: There are not enough draws computes to compute deciles. Skipping their computation.\n')
end

if np
    type = 'parameters';
    if TeX
        fid = TeXBegin(latexFolder, M_.fname, 1, type);
    end
    skipline()
    disp(type)
    disp(tit2)
    ip = nvx+nvn+ncx+ncn+1;
    for i=1:np
        if options_.mh_replic || (options_.load_mh_file && ~options_.load_results_after_load_mh) || ishssmc(options_)
            draws = getalldraws(ip);
            [post_mean, post_median, post_var, hpd_interval, post_deciles, density] = posterior_moments(draws, 1, options_.mh_conf_sig);
            name = bayestopt_.name{ip};
            oo_ = Filloo(oo_, name, type, post_mean, hpd_interval, post_median, post_var, post_deciles, density);
        else
            try
                name = bayestopt_.name{ip};
                [post_mean, hpd_interval, post_var] = Extractoo(oo_, name, type);
            catch
                draws = getalldraws(ip);
                [post_mean, post_median, post_var, hpd_interval, post_deciles, density] = posterior_moments(draws, 1, options_.mh_conf_sig);
                name = bayestopt_.name{ip};
                oo_ = Filloo(oo_, name, type, post_mean, hpd_interval, post_median, post_var, post_deciles, density);
            end
        end
        dprintf(pformat, header_width, name, bayestopt_.p1(ip),...
                post_mean, ...
                hpd_interval, ...
                pnames{bayestopt_.pshape(ip)+1}, ...
                bayestopt_.p2(ip));
        if TeX
            k = estim_params_.param_vals(i,1);
            name = M_.param_names_tex{k};
            TeXCore(fid, name, pnames{bayestopt_.pshape(ip)+1}, bayestopt_.p1(ip), bayestopt_.p2(ip), post_mean, sqrt(post_var), hpd_interval);
        end
        ip = ip+1;
    end
    if TeX
        TeXEnd(fid, 1, type);
    end
end

if nvx
    type = 'shocks_std';
    if TeX
        fid = TeXBegin(latexFolder, FileName,2, 'standard deviation of structural shocks');
    end
    ip = 1;
    skipline()
    disp('standard deviation of shocks')
    disp(tit2)
    for i=1:nvx
        if options_.mh_replic || (options_.load_mh_file && ~options_.load_results_after_load_mh)
            draws = getalldraws(ip);
            [post_mean, post_median, post_var, hpd_interval, post_deciles, density] = posterior_moments(draws, 1, options_.mh_conf_sig);
            k = estim_params_.var_exo(i,1);
            name = M_.exo_names{k};
            oo_ = Filloo(oo_, name, type, post_mean, hpd_interval, post_median, post_var, post_deciles, density);
            M_.Sigma_e(k,k) = post_mean*post_mean;
        else
            try
                k = estim_params_.var_exo(i,1);
                name = M_.exo_names{k};
                [post_mean, hpd_interval, post_var] = Extractoo(oo_, name, type);
            catch
                draws = getalldraws(ip);
                [post_mean, post_median, post_var, hpd_interval, post_deciles, density] = posterior_moments(draws, 1, options_.mh_conf_sig);
                k = estim_params_.var_exo(i,1);
                name = M_.exo_names{k};
                oo_ = Filloo(oo_, name, type, post_mean, hpd_interval, post_median, post_var, post_deciles, density);
                M_.Sigma_e(k,k) = post_mean*post_mean;
            end
        end
        dprintf(pformat, header_width, name, bayestopt_.p1(ip), post_mean, hpd_interval, pnames{bayestopt_.pshape(ip)+1}, bayestopt_.p2(ip));
        if TeX
            name = M_.exo_names_tex{k};
            TeXCore(fid,name, pnames{bayestopt_.pshape(ip)+1}, bayestopt_.p1(ip), bayestopt_.p2(ip), post_mean, sqrt(post_var), hpd_interval);
        end
        ip = ip+1;
    end
    if TeX
        TeXEnd(fid, 2, 'standard deviation of structural shocks');
    end
end

if nvn
    type = 'measurement_errors_std';
    if TeX
        fid = TeXBegin(latexFolder, FileName, 3, 'standard deviation of measurement errors');
    end
    skipline()
    disp('standard deviation of measurement errors')
    disp(tit2)
    ip = nvx+1;
    for i=1:nvn
        if options_.mh_replic || (options_.load_mh_file && ~options_.load_results_after_load_mh)
            draws = getalldraws(ip);
            [post_mean, post_median, post_var, hpd_interval, post_deciles, density] = posterior_moments(draws, 1, options_.mh_conf_sig);
            name = options_.varobs{estim_params_.nvn_observable_correspondence(i,1)};
            oo_ = Filloo(oo_, name, type, post_mean, hpd_interval, post_median, post_var, post_deciles, density);
        else
            try
                name = options_.varobs{estim_params_.nvn_observable_correspondence(i,1)};
                [post_mean,hpd_interval,post_var] = Extractoo(oo_,name,type);
            catch
                draws = getalldraws(ip);
                [post_mean, post_median, post_var, hpd_interval, post_deciles, density] = posterior_moments(draws,1,options_.mh_conf_sig);
                name = options_.varobs{estim_params_.nvn_observable_correspondence(i,1)};
                oo_ = Filloo(oo_,name,type,post_mean,hpd_interval,post_median,post_var,post_deciles,density);
            end
        end
        dprintf(pformat, header_width, name,bayestopt_.p1(ip), post_mean, hpd_interval, pnames{bayestopt_.pshape(ip)+1}, bayestopt_.p2(ip));
        if TeX
            k = estim_params_.var_endo(i,1);
            name = M_.endo_names_tex{k};
            TeXCore(fid, name, pnames{bayestopt_.pshape(ip)+1}, bayestopt_.p1(ip), bayestopt_.p2(ip), post_mean, sqrt(post_var), hpd_interval);
        end
        ip = ip+1;
    end
    if TeX
        TeXEnd(fid, 3, 'standard deviation of measurement errors');
    end
end

if ncx
    type = 'shocks_corr';
    if TeX
        fid = TeXBegin(latexFolder,FileName,4,'correlation of structural shocks');
    end
    skipline()
    disp('correlation of shocks')
    disp(tit2)
    ip = nvx+nvn+1;
    for i=1:ncx
        if options_.mh_replic || (options_.load_mh_file && ~options_.load_results_after_load_mh)
            draws = getalldraws(ip);
            [post_mean, post_median, post_var, hpd_interval, post_deciles, density] = posterior_moments(draws,1,options_.mh_conf_sig);
            k1 = estim_params_.corrx(i,1);
            k2 = estim_params_.corrx(i,2);
            name = sprintf('%s,%s', M_.exo_names{k1}, M_.exo_names{k2});
            NAME = sprintf('%s_%s', M_.exo_names{k1}, M_.exo_names{k2});
            oo_ = Filloo(oo_, NAME, type, post_mean, hpd_interval, post_median, post_var, post_deciles, density);
            M_.Sigma_e(k1,k2) = post_mean*sqrt(M_.Sigma_e(k1,k1)*M_.Sigma_e(k2,k2));
            M_.Sigma_e(k2,k1) = M_.Sigma_e(k1,k2);
        else
            try
                k1 = estim_params_.corrx(i,1);
                k2 = estim_params_.corrx(i,2);
                name = sprintf('%s,%s', M_.exo_names{k1}, M_.exo_names{k2});
                NAME = sprintf('%s_%s', M_.exo_names{k1}, M_.exo_names{k2});
                [post_mean,hpd_interval,post_var] = Extractoo(oo_, NAME, type);
            catch
                draws = getalldraws(ip);
                [post_mean, post_median, post_var, hpd_interval, post_deciles, density] = posterior_moments(draws, 1, options_.mh_conf_sig);
                k1 = estim_params_.corrx(i,1);
                k2 = estim_params_.corrx(i,2);
                name = sprintf('%s,%s', M_.exo_names{k1}, M_.exo_names{k2});
                NAME = sprintf('%s_%s', M_.exo_names{k1}, M_.exo_names{k2});
                oo_ = Filloo(oo_, NAME, type, post_mean, hpd_interval, post_median, post_var, post_deciles, density);
                M_.Sigma_e(k1,k2) = post_mean*sqrt(M_.Sigma_e(k1,k1)*M_.Sigma_e(k2,k2));
                M_.Sigma_e(k2,k1) = M_.Sigma_e(k1,k2);
            end
        end
        dprintf(pformat, header_width,name, bayestopt_.p1(ip), post_mean, hpd_interval, pnames{bayestopt_.pshape(ip)+1}, bayestopt_.p2(ip));
        if TeX
            name = sprintf('(%s,%s)', M_.exo_names_tex{k1}, M_.exo_names_tex{k2});
            TeXCore(fid, name, pnames{bayestopt_.pshape(ip)+1}, bayestopt_.p1(ip), bayestopt_.p2(ip), post_mean, sqrt(post_var), hpd_interval);
        end
        ip = ip+1;
    end
    if TeX
        TeXEnd(fid, 4, 'correlation of structural shocks');
    end
end

if ncn
    type = 'measurement_errors_corr';
    if TeX
        fid = TeXBegin(latexFolder, FileName, 5, 'correlation of measurement errors');
    end
    skipline()
    disp('correlation of measurement errors')
    disp(tit2)
    ip = nvx+nvn+ncx+1;
    for i=1:ncn
        if options_.mh_replic || (options_.load_mh_file && ~options_.load_results_after_load_mh)
            draws = getalldraws(ip);
            [post_mean, post_median, post_var, hpd_interval, post_deciles, density] = posterior_moments(draws, 1, options_.mh_conf_sig);
            k1 = estim_params_.corrn(i,1);
            k2 = estim_params_.corrn(i,2);
            name = sprintf('%s,%s', M_.endo_names{k1}, M_.endo_names{k2});
            NAME = sprintf('%s_%s', M_.endo_names{k1}, M_.endo_names{k2});
            oo_ = Filloo(oo_, NAME, type, post_mean, hpd_interval, post_median, post_var, post_deciles,density);
        else
            try
                k1 = estim_params_.corrn(i,1);
                k2 = estim_params_.corrn(i,2);
                name = sprintf('%s,%s', M_.endo_names{k1}, M_.endo_names{k2});
                NAME = sprintf('%s_%s', M_.endo_names{k1}, M_.endo_names{k2});
                [post_mean,hpd_interval,post_var] = Extractoo(oo_, NAME, type);
            catch
                draws = getalldraws(ip);
                [post_mean, post_median, post_var, hpd_interval, post_deciles, density] = posterior_moments(draws, 1, options_.mh_conf_sig);
                k1 = estim_params_.corrn(i,1);
                k2 = estim_params_.corrn(i,2);
                name = sprintf('%s,%s', M_.endo_names{k1}, M_.endo_names{k2});
                NAME = sprintf('%s_%s', M_.endo_names{k1}, M_.endo_names{k2});
                oo_ = Filloo(oo_, NAME, type, post_mean, hpd_interval, post_median, post_var, post_deciles, density);
            end
        end
        dprintf(pformat, header_width, name, bayestopt_.p1(ip), post_mean, hpd_interval, pnames{bayestopt_.pshape(ip)+1}, bayestopt_.p2(ip));
        if TeX
            name = sprintf('(%s,%s)', M_.endo_names_tex{k1}, M_.endo_names_tex{k2});
            TeXCore(fid, name, pnames{bayestopt_.pshape(ip)+1}, bayestopt_.p1(ip), bayestopt_.p2(ip), post_mean, sqrt(post_var), hpd_interval);
        end
        ip = ip+1;
    end
    if TeX
        TeXEnd(fid, 5, 'correlation of measurement errors');
    end
end


%
%% subfunctions:
%
function fid = TeXBegin(directory, fname, fnum, title)
TeXfile = [directory '/' fname '_Posterior_Mean_' int2str(fnum) '.tex'];
fidTeX = fopen(TeXfile,'w');
fprintf(fidTeX, '%% TeX-table generated by Dynare.\n');
fprintf(fidTeX, ['%% RESULTS FROM METROPOLIS HASTINGS (' title ')\n']);
fprintf(fidTeX, ['%% ' datestr(now, 0)]);
fprintf(fidTeX, ' \n');
fprintf(fidTeX, ' \n');
fprintf(fidTeX, '\\begin{center}\n');
fprintf(fidTeX, '\\begin{longtable}{llcccccc} \n');
fprintf(fidTeX, ['\\caption{Results from Metropolis-Hastings (' title ')}\n ']);
fprintf(fidTeX, ['\\label{Table:MHPosterior:' int2str(fnum)  '}\\\\\n']);
fprintf(fidTeX, '\\toprule \n');
fprintf(fidTeX, '  & \\multicolumn{3}{c}{Prior}  &  \\multicolumn{4}{c}{Posterior} \\\\\n');
fprintf(fidTeX, '  \\cmidrule(r{.75em}){2-4} \\cmidrule(r{.75em}){5-8}\n');
fprintf(fidTeX, '  & Dist. & Mean  & Stdev. & Mean & Stdev. & HPD inf & HPD sup\\\\\n');
fprintf(fidTeX, '\\midrule \\endfirsthead \n');
fprintf(fidTeX, '\\caption{(continued)}\\\\');
fprintf(fidTeX, '\\toprule \n');
fprintf(fidTeX, '  & \\multicolumn{3}{c}{Prior}  &  \\multicolumn{4}{c}{Posterior} \\\\\n');
fprintf(fidTeX, '  \\cmidrule(r{.75em}){2-4} \\cmidrule(r{.75em}){5-8}\n');
fprintf(fidTeX, '  & Dist. & Mean  & Stdev. & Mean & Stdev. & HPD inf & HPD sup\\\\\n');
fprintf(fidTeX, '\\midrule \\endhead \n');
fprintf(fidTeX, '\\bottomrule \\multicolumn{8}{r}{(Continued on next page)} \\endfoot \n');
fprintf(fidTeX, '\\bottomrule \\endlastfoot \n');
fid = fidTeX;


function TeXCore(fid, name, shape, priormean, priorstd, postmean, poststd, hpd)
fprintf(fid,'$%s$ & %s & %7.3f & %6.4f & %7.3f& %6.4f & %7.4f & %7.4f \\\\ \n',...
        name, ...
        shape, ...
        priormean, ...
        priorstd, ...
        postmean, ...
        poststd, ...
        hpd(1), ...
        hpd(2));


function TeXEnd(fid, fnum, title)
fprintf(fid, '\\end{longtable}\n ');
fprintf(fid, '\\end{center}\n');
fprintf(fid, '%% End of TeX file.\n');
fclose(fid);


function oo = Filloo(oo, name, type, postmean, hpdinterval, postmedian, postvar, postdecile, density)
oo.posterior_mean.(type).(name) = postmean;
oo.posterior_hpdinf.(type).(name) = hpdinterval(1);
oo.posterior_hpdsup.(type).(name) = hpdinterval(2);
oo.posterior_median.(type).(name) = postmedian;
oo.posterior_variance.(type).(name) = postvar;
oo.posterior_std.(type).(name) = sqrt(postvar);
oo.posterior_deciles.(type).(name) = postdecile;
oo.posterior_density.(type).(name) = density;

function [post_mean,hpd_interval,post_var] = Extractoo(oo, name, type)
hpd_interval = zeros(2,1);
post_mean = oo.posterior_mean.(type).(name);
hpd_interval(1) = oo.posterior_hpdinf.(type).(name);
hpd_interval(2) = oo.posterior_hpdsup.(type).(name);
post_var = oo.posterior_variance.(type).(name);
