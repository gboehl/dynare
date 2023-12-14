function indmcf = monte_carlo_filtering_analysis(lpmat, ibeha, inobeha, options_mcf, M_, options_, bayestopt_, estim_params_)
% indmcf = monte_carlo_filtering_analysis(lpmat, ibeha, inobeha, options_mcf, M_, options_, bayestopt_, estim_params_)
% Inputs:
% - lpmat               [double]        Monte Carlo matrix
% - ibeha               [integer]       index of behavioural runs
% - inobeha             [integer]       index of non-behavioural runs
% - options_gsa_        [structure]     GSA options_
% - M_                  [structure]     describing the model
% - options_            [structure]     describing the options
% - bayestopt_          [structure]     describing the priors
% - estim_params_       [structure]     characterizing parameters to be estimated
%
% Outputs:
% - indmcf              [double]        results of matrix

% Written by Marco Ratto
% Joint Research Centre, The European Commission,
% marco.ratto@ec.europa.eu
%

% Copyright © 2014 European Commission
% Copyright © 2016-2023 Dynare Team
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

pvalue_ks = options_mcf.pvalue_ks;
pvalue_corr = options_mcf.pvalue_corr;
alpha2 = options_mcf.alpha2;
param_names = options_mcf.param_names;

if options_.TeX
    if ~isfield(options_mcf,'param_names_tex')
        param_names_tex = options_mcf.param_names;
    else
        param_names_tex = options_mcf.param_names_tex;
    end
else
    param_names_tex = strrep(options_mcf.param_names,'_','\_');
end
amcf_name = options_mcf.amcf_name;
amcf_title = options_mcf.amcf_title;
beha_title = options_mcf.beha_title;
nobeha_title = options_mcf.nobeha_title;
if options_.TeX
    beha_title_latex = options_mcf.beha_title_latex;
    nobeha_title_latex = options_mcf.nobeha_title_latex;
end
title = options_mcf.title;
fname_ = options_mcf.fname_;
xparam1=[];
if isfield(options_mcf,'xparam1')
    xparam1=options_mcf.xparam1;
end
OutputDirectoryName = options_mcf.OutputDirectoryName;

[proba, dproba] = gsa.stability_mapping_univariate(lpmat, ibeha, inobeha, [],fname_, options_, bayestopt_.name, estim_params_,0);
indmcf=find(proba<pvalue_ks);
[~,jtmp] = sort(proba(indmcf),1,'ascend');
indmcf = indmcf(jtmp);
if ~isempty(indmcf)
    skipline()
    headers = {'Parameter','d-stat','p-value'};
    labels = param_names(indmcf);
    data_mat=[dproba(indmcf) proba(indmcf)];
    options_temp.noprint=0;
    dyntable(options_temp,['Smirnov statistics in driving ', title],headers,labels,data_mat,size(labels,2)+2,16,3);
    if options_.TeX
        labels_TeX=param_names_tex(indmcf);
        M_temp.dname=OutputDirectoryName ;
        M_temp.fname=fname_;
        dyn_latex_table(M_temp,options_temp,['Smirnov statistics in driving ', strrep(title,'_','\\_')],amcf_name,headers,labels_TeX,data_mat,size(labels,2)+2,16,6);
    end
end

if length(ibeha)>10 && length(inobeha)>10
    if options_.TeX
        indcorr1 = gsa.stability_mapping_bivariate(lpmat(ibeha,:),alpha2, pvalue_corr, M_, options_, bayestopt_, estim_params_, beha_title, beha_title_latex);
        indcorr2 = gsa.stability_mapping_bivariate(lpmat(inobeha,:),alpha2, pvalue_corr, M_, options_, bayestopt_, estim_params_, nobeha_title, nobeha_title_latex);
    else
        indcorr1 = gsa.stability_mapping_bivariate(lpmat(ibeha,:),alpha2, pvalue_corr, M_, options_, bayestopt_, estim_params_, beha_title);
        indcorr2 = gsa.stability_mapping_bivariate(lpmat(inobeha,:),alpha2, pvalue_corr, M_, options_, bayestopt_, estim_params_, nobeha_title);
    end    
    indcorr = union(indcorr1(:), indcorr2(:));
    indcorr = indcorr(~ismember(indcorr(:),indmcf));
    indmcf = [indmcf(:); indcorr(:)];
end
if ~isempty(indmcf) && ~options_.nograph
    skipline()
    xx=[];
    if ~ isempty(xparam1)
        xx=xparam1(indmcf); 
    end
    if options_.TeX
        gsa.scatter_mcf(lpmat(ibeha,indmcf),lpmat(inobeha,indmcf), param_names_tex(indmcf), ...
            '.', [fname_,'_',amcf_name], OutputDirectoryName, amcf_title,xx, options_, ...
            beha_title, nobeha_title, beha_title_latex, nobeha_title_latex)
    else
        gsa.scatter_mcf(lpmat(ibeha,indmcf),lpmat(inobeha,indmcf), param_names_tex(indmcf), ...
            '.', [fname_,'_',amcf_name], OutputDirectoryName, amcf_title,xx, options_, ...
            beha_title, nobeha_title)
    end
end
