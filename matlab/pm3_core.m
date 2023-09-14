function myoutput=pm3_core(myinputs,fpar,nvar,whoiam, ThisMatlab)
% myoutput=pm3_core(myinputs,fpar,nvar,whoiam, ThisMatlab)
% PARALLEL CONTEXT
% Core functionality for pm3.m function, which can be parallelized.

% INPUTS
%   o myimput            [struc]     The mandatory variables for local/remote
%                                    parallel computing obtained from prior_posterior_statistics.m
%                                    function.
%   o fpar and nvar      [integer]   first variable and number of variables
%   o whoiam             [integer]   In concurrent programming a modality to refer to the different threads running in parallel is needed.
%                                    The integer whoaim is the integer that
%                                    allows us to distinguish between them. Then it is the index number of this CPU among all CPUs in the
%                                    cluster.
%   o ThisMatlab         [integer]   Allows us to distinguish between the
%                                    'main' Matlab, the slave Matlab worker, local Matlab, remote Matlab,
%                                     ... Then it is the index number of this slave machine in the cluster.
%
% OUTPUTS
% o myoutput              [struct]  Contains file names
%
%
% SPECIAL REQUIREMENTS.
%   None.

% Copyright © 2007-2023 Dynare Team
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

if nargin<4
    whoiam=0;
end

% Reshape 'myinputs' for local computation.
% In order to avoid confusion in the name space, the instruction struct2local(myinputs) is replaced by:

tit1=myinputs.tit1;
nn=myinputs.nn;
n2=myinputs.n2;
Distrib=myinputs.Distrib;
varlist=myinputs.varlist;

MaxNumberOfPlotsPerFigure=myinputs.MaxNumberOfPlotsPerFigure;
name3=myinputs.name3;
tit2=myinputs.tit2;
Mean=myinputs.Mean;
options_.TeX=myinputs.TeX;
options_.nodisplay=myinputs.nodisplay;
options_.graph_format=myinputs.graph_format;
fname=myinputs.fname;
dname=myinputs.dname;

if whoiam
    Parallel=myinputs.Parallel;
end

if options_.TeX
    varlist_TeX=myinputs.varlist_TeX;
end

if whoiam
    prct0={0,whoiam,Parallel(ThisMatlab)};
    h = dyn_waitbar(prct0,'Parallel plots pm3 ...');
end

figunumber = 0;
subplotnum = 0;
hh_fig = dyn_figure(options_.nodisplay,'Name',[tit1 ' ' int2str(figunumber+1)]);
RemoteFlag = 0;
if whoiam
    if Parallel(ThisMatlab).Local ==0
        RemoteFlag=1;
    end
end

OutputFileName = {};

for i=fpar:nvar
    if max(abs(Mean(:,i))) > 10^(-6)
        subplotnum = subplotnum+1;
        set(0,'CurrentFigure',hh_fig);
        subplot(nn,nn,subplotnum);
        if ~(all(all(Distrib(:,:,i)))>0 || ~all(all(Distrib(:,:,i)))<0)
            plot([1 n2],[0 0],'-r','linewidth',0.5);
        end
        hold on
        for k = 1:9
            plot(1:n2,squeeze(Distrib(k,:,i)),'-g','linewidth',0.5);
        end
        plot(1:n2,Mean(:,i),'-k','linewidth',1);
        xlim([1 n2]);
        hold off;
        if options_.TeX
            title(['$' varlist_TeX{i,:} '$'],'Interpreter','latex')            
        else
            title(varlist(i,:),'Interpreter','none')
        end       
        yticklabels=get(gca,'yticklabel');
        if size(char(yticklabels),2)>5 %make sure yticks do not screw up figure
            yticks=get(gca,'ytick');
            for ii=1:length(yticks)
                yticklabels_new{ii,1}=sprintf('%4.3f',yticks(ii));
            end
            set(gca,'yticklabel',yticklabels_new)
        end
    end

    if whoiam
        if Parallel(ThisMatlab).Local==0
            DirectoryName = CheckPath('Output',dname);
        end
    end

    if subplotnum == MaxNumberOfPlotsPerFigure || i == nvar
        dyn_saveas(hh_fig,[dname '/Output/'  fname '_' name3 '_' tit2{i}],options_.nodisplay,options_.graph_format);
        if RemoteFlag==1
            OutputFileName = [OutputFileName; {[dname, filesep, 'Output',filesep], [fname '_' name3 '_' deblank(tit2(i,:)) '.*']}];
        end
        subplotnum = 0;
        figunumber = figunumber+1;
        if (i ~= nvar)
            hh_fig = dyn_figure(options_.nodisplay,'Name',[name3 ' ' int2str(figunumber+1)]);
        end
    end

    if whoiam
        dyn_waitbar((i-fpar+1)/(nvar-fpar+1),h);
    end
end

if whoiam
    dyn_waitbar_close(h);
end
myoutput.OutputFileName=OutputFileName;
