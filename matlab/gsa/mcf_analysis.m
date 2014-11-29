function mcf_analysis(lpmat, ibeha, inobeha, options_mcf, DynareOptions)
%
% Written by Marco Ratto
% Joint Research Centre, The European Commission,
% marco.ratto@jrc.ec.europa.eu
%

% Copyright (C) 2014 European Commission
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
% along with Dynare.  If not, see <http://www.gnu.org/licenses/>.

pvalue_ks = options_mcf.pvalue_ks;
pvalue_corr = options_mcf.pvalue_corr;
alpha2 = options_mcf.alpha2;
param_names = options_mcf.param_names;
amcf_name = options_mcf.amcf_name;
amcf_title = options_mcf.amcf_title;
beha_title = options_mcf.beha_title;
nobeha_title = options_mcf.nobeha_title;
title = options_mcf.title;
fname_ = options_mcf.fname_;
OutputDirectoryName = options_mcf.OutputDirectoryName;

[proba, dproba] = stab_map_1(lpmat, ibeha, inobeha, [],0);
%         indindet=find(dproba>ksstat);
indmcf=find(proba<pvalue_ks);
[tmp,jtmp] = sort(proba(indmcf),2,'ascend');
indmcf = indmcf(jtmp);
if ~isempty(indmcf)
    disp(['Smirnov statistics in driving ', title])
    for j=1:length(indmcf),
        disp([param_names(indmcf(j),:),'   d-stat = ', num2str(dproba(indmcf(j)),'%1.3f'),'   p-value = ', num2str(proba(indmcf(j)),'%1.3f')])
    end
    skipline()
end
if length(ibeha)>10 && length(inobeha)>10,
    indcorr1 = stab_map_2(lpmat(ibeha,:),alpha2, pvalue_corr, beha_title);
    indcorr2 = stab_map_2(lpmat(inobeha,:),alpha2, pvalue_corr, nobeha_title);
    indcorr = union(indcorr1(:), indcorr2(:));
    indcorr = indcorr(~ismember(indcorr(:),indmcf));
    indmcf = [indmcf(:); indcorr(:)];
end
if ~isempty(indmcf)
    skipline()
    scatter_mcf(lpmat(ibeha,indmcf),lpmat(inobeha,indmcf), param_names(indmcf,:), ...
        '.', [fname_,'_',amcf_name], OutputDirectoryName, amcf_title,[], DynareOptions, ...
        beha_title, nobeha_title)
end
