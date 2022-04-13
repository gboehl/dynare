%function msstart2(options_)
%   This .m file is called for point graphics or error bands and
%   It starts for both data and Bayesian estimation (when IxEstimat==0,
%          no estimation but only data analysis), which allows you to set
% (a) available data range,
% (b) sample range,
% (c) rearrangement of actual data such as mlog, qg, yg
% (d) conditions of shocks 'Cms',
% (c) conditions of variables 'nconstr',
% (e) soft conditions 'nbancon,'
% (f) produce point conditional forecast (at least conditions on variables).
%
% February 2004

% Copyright © 2011-2017 Dynare Team
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

% ** ONLY UNDER UNIX SYSTEM
%path(path,'/usr2/f1taz14/mymatlab')

%addpath('C:\SoftWDisk\MATLAB6p5\toolbox\cstz')


msstart_setup
%===========================================
% Exordium II
%===========================================
%options_.ms.ncsk = 0;   % conditional directly on shoocks.  Unlike Cms, not on variables that back out shocks
%options_.ms.nstd = 6;   % number of standard deviations to cover the range in which distributions are put to bins
%options_.ms.ninv = 1000;   % the number of bins for grouping, say, impulse responses
%options_.ms.Indxcol = [1:nvar];  % a vector of random columns in which MC draws are made.
%
%options_.ms.indxparr = 1;  % 1: parameters random;  0: no randomness in parameters
% Note, when 0, there is no effect from the values of options_.ms.IndxAp, options_.ms.Aband, etc.
%options_.ms.indxovr = 0;   % 1: distributions for other variables of interest; 0: no distribution.
% Example:  joint distribution of a(1) and a(2).  Only for specific purposes
%options_.ms.Aband = 1;     % 1: error bands with only A0 and A+ random.
%options_.ms.IndxAp = 1;    % 1: generate draws of A+; 0: no such draws.
% Note: when options_.ms.IndxAp=0, there is no effect from the values of options_.ms.options_.ms.options_.ms.options_.ms.indximf, IndxFore,
%         or options_.ms.apband.
%options_.ms.apband = 1;    % 1: error bands for A+; 0: no error bands for A+.
%*** The following (impulse responses and forecasts) is used only if options_.ms.IndxAp=1
%options_.ms.indximf = 1;   % 1: generate draws of impulse responses; 0: no such draws (thus no effect
%      from options_.ms.imfband)
%options_.ms.imfband = 1;   % 1: error bands for impulse responses; 0: no error bands
%options_.ms.indxfore = 0;  % 1: generate draws of forecasts; 0: no such draws (thus no effect from options_.ms.foreband)
%options_.ms.foreband = 0;  % 1: error bands for out-of-sample forecasts; 0: no error bands
%
%options_.ms.indxgforhat = 1;  % 1: plot unconditoinal forecasts; 0: no such plot
rnum = nvar;      % number of rows in the graph
cnum = 1;      % number of columns in the graph
if rnum*cnum<nvar
    warning('rnum*cnum must be at least as large as nvar')
    disp('Hit any key to continue, or ctrl-c to abort')
    pause
end
%options_.ms.indxgimfhat = 1;  % 1: plot ML impulse responses; 0: no plot
%options_.ms.indxestima = 1;  %1: ML estimation; 0: no estimation and data only
%
IndxNmlr = [1 0 0 0 0 0];  % imported by nmlzvar.m
                           % Index for which normalization rule to choose
                           % Only one of the elments in IndxNmlr can be non-zero
                           % IndxNmlr(1): ML A distance rule (supposed to be the best)
                           % IndxNmlr(2): ML Ahat distance rule (to approximate IndxNmlr(1))
                           % IndxNmlr(3): ML Euclidean distance rule (not invariant to scale)
                           % IndxNmlr(4): Positive diagonal rule
                           % IndxNmlr(5): Positive inv(A) diagonal rule (if ~IndxNmlr(5), no need for A0inu,
                           %                                      so let A0inu=[])
                           % IndxNmlr(6): Assigned postive rule (such as off-diagonal elements).  Added 1/3/00


%%%%----------------------------------------
% Hard conditions directly on variables
%
%options_.ms.indxgdls = 1;  % 1: graph point forecast on variables; 0: disable
nconstr1=nfqm;      % number of the 1st set of constraints
nconstr2=options_.forecast ;     % number of the 2nd set of constraints
nconstr=nconstr1+nconstr2;   % q: 4 years -- 4*12 months.
                             % When 0, no conditions directly on variables <<>>
nconstr=0 ;  %6*nconstr1;
options_.ms.eq_ms = [];      % location of MS equation; if [], all shocks
PorR = [4*ones(nconstr1,1);2*ones(nconstr1,1);3*ones(nconstr1,1)];   % the variable conditioned.  1: Pcm; 3: FFR; 4: CPI
PorR = [PorR;1*ones(nconstr1,1);5*ones(nconstr1,1);6*ones(nconstr1,1)];
%PorR = [3 5];   % the variable conditioned.  3: FFR; 5: CPI
%PorR = 3;

%%%%----------------------------------------
% Conditions directly on future shocks
%
%options_.ms.cms = 0     % 1: condition on ms shocks; 0: disable this and "fidcnderr.m" gives
%   unconditional forecasts if nconstr = 0 as well;  <<>>
%options_.ms.ncms = 0;   % number of the stance of policy; 0 if no tightening or loosening
%options_.ms.eq_cms = 1;  % location of MS shocks
options_.ms.tlindx = 1*ones(1,options_.ms.ncms);  % 1-by-options_.ms.ncms vector; 1: tightening; 0: loosen
options_.ms.tlnumber = [0.5 0.5 0 0];  %94:4 % [2 2 1.5 1.5]; %79:9  %[1.5 1.5 1 1]; 90:9
                                       % 1-by-options_.ms.ncms vector; cut-off point for MS shocks
TLmean = zeros(1,options_.ms.ncms);
% unconditional, i.e., 0 mean, for the final report in the paper
if options_.ms.cms
    options_.ms.eq_ms = [];
    % At least at this point, it makes no sense to have DLS type of options_.ms.eq_ms; 10/12/98
    if all(isfinite(options_.ms.tlnumber))
        for k=1:options_.ms.ncms
            TLmean(k) = lcnmean(options_.ms.tlnumber(k),options_.ms.tlindx(k));
            % shock mean magnitude. 1: tight; 0: loose
            % Never used for any subsequent computation but
            %   simply used for the final report in the paper.
            %options_.ms.tlnumber(k) = fzero('lcutoff',0,[],[],TLmean(k))
            % get an idea about the cutoff point given TLmean instead

        end
    end
else
    options_.ms.ncms = 0;   % only for the use of the graph by msprobg.m
    options_.ms.tlnumber = NaN*ones(1,options_.ms.ncms);
    % -infinity, only for the use of the graph by msprobg.m
end


%%%%----------------------------------------
% Soft conditions on variables
%
%cnum = 0  % # of band condtions; when 0, disable this option
% Note (different from "fidencon") that each condition corres. to variable
%options_.ms.banact = 1;    % 1: use infor on actual; 0:  preset without infor on actual
if cnum
    banindx = cell(cnum,1);  % index for each variable or conditon
    banstp = cell(cnum,1);    % steps:  annual in general
    banvar = zeros(cnum,1);    % varables:  annual in general
    banval = cell(cnum,1);    % band value (each variable occupy a cell)
    badval{1} = zeros(length(banstp{1}),2);   % 2: lower or higher bound

    banstp{1} = 1:4;      % 3 or 4 years
    banvar(1) = 3;      % 3: FFR; 5: CPI
    if ~options_.ms.banact
        for i=1:length(banstp{1})
            banval{1}(i,:) = [5.0 10.0];
        end
    end
end
%
pause(1)
skipline()
disp('For uncondtional forecasts, set nconstr = options_.ms.cms = cnum = 0')
pause(1)
%
%=================================================
%====== End of exordium ==========================
%=================================================





%(1)--------------------------------------
% Further data analysis
%(1)--------------------------------------
%
if (options_.ms.freq==12)
    nStart=(yrStart-options_.ms.initial_year )*12+qmStart-options_.ms.initial_subperiod ;  % positive number of months at the start
    nEnd=(yrEnd-options_.ms.final_year )*12+qmEnd-options_.ms.final_subperiod ;     % negative number of months towards end
elseif (options_.ms.freq==4)
    nStart=(yrStart-options_.ms.initial_year )*4+qmStart-options_.ms.initial_subperiod ;  % positive number of months at the start
    nEnd=(yrEnd-options_.ms.final_year )*4+qmEnd-options_.ms.final_subperiod ;     % negative number of months towards end
elseif (options_.ms.freq==1)
    nStart=(yrStart-options_.ms.initial_year )*1+qmStart-options_.ms.initial_subperiod ;  % positive number of months at the start
    nEnd=(yrEnd-options_.ms.final_year )*1+qmEnd-options_.ms.final_subperiod ;     % negative number of months towards end
else
    error('Error: this code is only good for monthly/quarterly/yearly data!!!')
    return
end
%
if nEnd>0 || nStart<0
    disp('Warning: this particular sample consider is out of bounds of the data!!!')
    return
end
%***  Note, both xdgel and xdata have the same start with the specific sample
xdgel=options_.data(nStart+1:nData+nEnd,options_.ms.vlist);
% gel: general options_.data within sample (nSample)
if ~(nSample==size(xdgel,1))
    warning('The sample size (including options_.ms.nlags ) and data are incompatible')
    disp('Check to make sure nSample and size(xdgel,1) are the same')
    return
end
%
baddata = find(isnan(xdgel));
if ~isempty(baddata)
    warning('Some data for this selected sample are actually unavailable.')
    disp('Hit any key to continue, or ctrl-c to abort')
    pause
end
%
if options_.ms.initial_subperiod ==1
    yrB = options_.ms.initial_year ; qmB = options_.ms.initial_subperiod ;
else
    yrB = options_.ms.initial_year +1; qmB = 1;
end
yrF = options_.ms.final_year ; qmF = options_.ms.final_subperiod ;
[Mdate,tmp] = fn_calyrqm(options_.ms.freq,[options_.ms.initial_year  options_.ms.initial_subperiod ],[options_.ms.final_year options_.ms.final_subperiod ]);
xdatae=[Mdate options_.data(1:nData,options_.ms.vlist)];
% beyond sample into forecast horizon until the end of the data options_.ms.final_year :options_.ms.final_subperiod
% Note: may contain NaN data.  So must be careful about its use

%=========== Obtain prior-period, period-to-last period, and annual growth rates
[yactyrge,yactyre,yactqmyge,yactqmge,yactqme] = fn_datana(xdatae,options_.ms.freq,options_.ms.log_var,options_.ms.percent_var,[yrB qmB],[yrF qmF]);
qdates = zeros(size(yactqmyge,1),1);
for ki=1:length(qdates)
    qdates(ki) = yactqmyge(1,1) + (yactqmyge(1,2)+ki-2)/options_.ms.freq;
end
for ki=1:nvar
    figure
    plot(qdates, yactqmyge(:,2+ki)/100)
    xlabel(options_.ms.varlist{ki})
end
save outactqmygdata.prn yactqmyge -ascii



%===========  Write the output on the screen or to a file in an organized way ==============
%disp([sprintf('%4.0f %2.0f %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f\n',yactyrge')])
spstr1 = 'disp([sprintf(';
spstr2 = '%4.0f %2.0f';
yactyrget=yactyrge';
for ki=1:length(options_.ms.vlist)
    if ki==length(options_.ms.vlist)
        spstr2 = [spstr2 ' %8.3f\n'];
    else
        spstr2 = [spstr2 ' %8.3f'];
    end
end
spstr = [spstr1 'spstr2' ', yactyrget)])'];
eval(spstr)

%
fid = fopen('outyrqm.prn','w');
%fprintf(fid,'%4.0f %2.0f %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f\n',yactyrge');
fpstr1 = 'fprintf(fid,';
fpstr2 = '%4.0f %2.0f';
for ki=1:nvar
    if ki==nvar
        fpstr2 = [fpstr2 ' %8.3f\n'];
    else
        fpstr2 = [fpstr2 ' %8.3f'];
    end
end
fpstr = [fpstr1 'fpstr2' ', yactyrget);'];
eval(fpstr)
fclose(fid);



if options_.ms.indxestima
    %(2)----------------------------------------------------------------------------
    % Estimation
    % ML forecast and impulse responses
    % Hard or soft conditions for conditional forecasts
    %(2)----------------------------------------------------------------------------
    %
    %* Arranged data information, WITHOUT dummy obs when 0 after mu is used.  See fn_rnrprior_covres_dobs.m for using the dummy
    %    observations as part of an explicit prior.
    [xtx,xty,yty,fss,phi,y,ncoef,xr,Bh] = fn_dataxy(nvar,options_.ms.nlags ,xdgel,mu,0,nexo);
    if qmStart+options_.ms.nlags -options_.ms.dummy_obs >0
        qmStartEsti = rem(qmStart+options_.ms.nlags -options_.ms.dummy_obs ,options_.ms.freq);   % dummy observations are included in the sample.
        if (~qmStartEsti)
            qmStartEsti = options_.ms.freq;
        end
        yrStartEsti = yrStart + floor((qmStart+options_.ms.nlags -options_.ms.dummy_obs )/(options_.ms.freq+0.01));
        % + 0.01 (or any number < 1)  is used so that qmStart+options_.ms.nlags -options_.ms.dummy_obs ==?*options_.ms.freq doesn't give us an extra year forward.
    else
        qmStartEsti = options_.ms.freq + rem(qmStart+options_.ms.nlags -options_.ms.dummy_obs ,options_.ms.freq);   % dummy observations are included in the sample.
        if (qmStart+options_.ms.nlags -options_.ms.dummy_obs ==0)
            yrStartEsti = yrStart - 1;   % one year back.
        else
            yrStartEsti = yrStart + floor((qmStart+options_.ms.nlags -options_.ms.dummy_obs )/(options_.ms.freq-0.01));
            % - 0.01 (or any number < 1)  is used so that qmStart+options_.ms.nlags -options_.ms.dummy_obs ==-?*options_.ms.freq give us an extra year back.
        end
    end
    dateswd = fn_dataext([yrStartEsti qmStartEsti],[yrEnd qmEnd],xdatae(:,[1:2]));  % dates with dummies
    phie = [dateswd phi];
    ye = [dateswd y];

    %* Obtain linear restrictions
    [Uiconst,Viconst,n0,np,ixmC0Pres] = feval(options_.ms.restriction_fname,nvar,nexo,options_.ms );
    if min(n0)==0
        skipline()
        warning('A0: restrictions in dlrprior.m give no free parameter in one of equations')
        disp('Press ctrl-c to abort')
        pause
    elseif min(np)==0
        skipline()
        warning('Ap: Restrictions in dlrprior.m give no free parameter in one of equations')
        disp('Press ctrl-c to abort')
        pause
    end

    if options_.ms.contemp_reduced_form
        Uiconst=cell(nvar,1); Viconst=cell(ncoef,1);
        for kj=1:nvar
            Uiconst{kj} = eye(nvar);  Viconst{kj} = eye(ncoef);
        end
    end

    if options_.ms.bayesian_prior
        %*** Obtains asymmetric prior (with no linear restrictions) with dummy observations as part of an explicit prior (i.e,
        %      reflected in Hpmulti and Hpinvmulti).  See Forecast II, pp.69a-69b for details.
        if 1  % Liquidity effect prior on both MS and MD equations.
            [Pi,H0multi,Hpmulti,H0invmulti,Hpinvmulti] = fn_rnrprior_covres_dobs(nvar,options_.ms.freq,options_.ms.nlags ,xdgel,mu,indxDummy,hpmsmd,indxmsmdeqn);
        else
            [Pi,H0multi,Hpmulti,H0invmulti,Hpinvmulti] = fn_rnrprior(nvar,options_.ms.freq,options_.ms.nlags ,xdgel,mu);
        end

        %*** Combines asymmetric prior with linear restrictions
        [Ptld,H0invtld,Hpinvtld] = fn_rlrprior(Uiconst,Viconst,Pi,H0multi,Hpmulti,nvar);

        %*** Obtains the posterior matrices for estimation and inference
        [Pmat,H0inv,Hpinv] = fn_rlrpostr(xtx,xty,yty,Ptld,H0invtld,Hpinvtld,Uiconst,Viconst);

        if options_.ms.contemp_reduced_form
            %*** Obtain the ML estimate
            A0hatinv = chol(H0inv{1}/fss);   % upper triangular but lower triangular choleski
            A0hat=inv(A0hatinv);
            a0indx = find(A0hat);
        else
            %*** Obtain the ML estimate
            %   load idenml
            x = 10*rand(sum(n0),1);
            H0 = eye(sum(n0));
            crit = 1.0e-9;
            nit = 10000;
            %
            [fhat,xhat,grad,Hhat,itct,fcount,retcodehat] = ...
                csminwel('fn_a0freefun',x,H0,'fn_a0freegrad',crit,nit,Uiconst,nvar,n0,fss,H0inv);

            A0hat = fn_tran_b2a(xhat,Uiconst,nvar,n0)
            A0hatinv = inv(A0hat);
            fhat
            xhat
            grad
            itct
            fcount
            retcodehat
            save outm.mat xhat A0hat A0hatinv grad fhat itct itct fcount retcodehat
        end
    else
        %*** Obtain the posterior matrices for estimation and inference
        [Pmat,H0inv,Hpinv] = fn_dlrpostr(xtx,xty,yty,Uiconst,Viconst);

        if options_.ms.contemp_reduced_form
            %*** Obtain the ML estimate
            A0hatinv = chol(H0inv{1}/fss);   % upper triangular but lower triangular choleski
            A0hat=inv(A0hatinv);
            a0indx = find(A0hat);
        else
            %*** Obtain the ML estimate
            %   load idenml
            x = 10*rand(sum(n0),1);
            H0 = eye(sum(n0));
            crit = 1.0e-9;
            nit = 10000;
            %
            [fhat,xhat,grad,Hhat,itct,fcount,retcodehat] = ...
                csminwel('fn_a0freefun',x,H0,'fn_a0freegrad',crit,nit,Uiconst,nvar,n0,fss,H0inv);

            A0hat = fn_tran_b2a(xhat,Uiconst,nvar,n0)
            A0hatinv = inv(A0hat);
            fhat
            xhat
            grad
            itct
            fcount
            retcodehat
            save outm.mat xhat A0hat A0hatinv grad fhat itct itct fcount retcodehat
        end
    end

    %**** impulse responses
    swish = A0hatinv;       % each column corresponds to an equation
    if options_.ms.contemp_reduced_form
        xhat = A0hat(a0indx);
        Bhat=Pmat{1};
        Fhat = Bhat*A0hat
        ghat = NaN;
    else
        xhat = fn_tran_a2b(A0hat,Uiconst,nvar,n0);
        [Fhat,ghat] = fn_gfmean(xhat,Pmat,Viconst,nvar,ncoef,n0,np);
        if options_.ms.cross_restrictions
            Fhatur0P = Fhat;  % ur: unrestriced across A0 and A+
            for ki = 1:size(ixmC0Pres,1)   % loop through the number of equations in which
                                           % cross-A0-A+ restrictions occur. See St. Louis Note p.5.
                ixeq = ixmC0Pres{ki}(1,1);   % index for the jth equation in consideration.
                Lit = Viconst{ixeq}(ixmC0Pres{ki}(:,2),:);  % transposed restriction matrix Li
                                                            % V_j(i,:) in f_j(i) = V_j(i,:)*g_j
                ci = ixmC0Pres{ki}(:,4) .* A0hat(ixmC0Pres{ki}(:,3),ixeq);
                % s * a_j(h) in the restriction f_j(i) = s * a_j(h).
                LtH = Lit/Hpinv{ixeq};
                HLV = LtH'/(LtH*Lit');
                gihat = Viconst{ixeq}'*Fhatur0P(:,ixeq);
                Fhat(:,ixeq) = Viconst{ixeq}*(gihat + HLV*(ci-Lit*gihat));
            end
        end
        Fhat
        Bhat = Fhat/A0hat;   % ncoef-by-nvar reduced form lagged parameters.
    end
    nn = [nvar options_.ms.nlags  imstp];
    imfhat = fn_impulse(Bhat,swish,nn);    % in the form that is congenial to RATS
    imf3hat=reshape(imfhat,size(imfhat,1),nvar,nvar);
    % imf3: row--steps, column--nvar responses, 3rd dimension--nvar shocks
    imf3shat=permute(imf3hat,[1 3 2]);
    % imf3s: permuted so that row--steps, column--nvar shocks,
    %                                3rd dimension--nvar responses
    % Note: reshape(imf3s(1,:,:),nvar,nvar) = A0in  (columns -- equations)
    if options_.ms.indxgimfhat
        figure
    end
    scaleout = fn_imcgraph(imfhat,nvar,imstp,xlab,ylab,options_.ms.indxgimfhat);
    imfstd = max(abs(scaleout)');   % row: nvar (largest number); used for standard deviations

    %
    %  %**** save stds. of both data and impulse responses in idfile1
    %  temp = [std(yactqmyge(:,3:end)); std(yactyrge(:,3:end)); imfstd];  %<<>>
    %  save idenyimstd.prn temp -ascii   % export forecast and impulse response to the file "idenyimstd.prn", 3-by-nvar
    %  %
    %  %**** save stds. of both data and impulse responses in idfile1
    %  temp = [std(yactqmyge(:,3:end)); std(yactyrge(:,3:end)); imfstd];  %<<>>
    %  save idenyimstd.prn temp -ascii   % export forecast and impulse response to the file "idenyimstd.prn", 3-by-nvar
    %  if options_.ms.indxparr
    %     idfile1='idenyimstd';
    %  end

    %=====================================
    % Now, out-of-sample forecasts. Note: Hm1t does not change with A0.
    %=====================================
    %
    % * updating the last row of X (phi) with the current (last row of) y.
    tcwx = nvar*options_.ms.nlags ;  % total coefficeint without exogenous variables
    phil = phi(size(phi,1),:);
    phil(nvar+1:tcwx) = phil(1:tcwx-nvar);
    phil(1:nvar) = y(end,:);
    %*** exogenous variables excluding constant terms
    if (nexo>1)
        Xexoe = fn_dataext([yrEnd qmEnd],[yrEnd qmEnd],xdatae(:,[1:2 2+nvar+1:2+nvar+nexo-1]));
        phil(1,tcwx+1:tcwx+nexo-1) = Xexoe(1,3:end);
    end
    %
    %*** ML unconditional point forecast
    nn = [nvar options_.ms.nlags  nfqm];
    if nexo<2
        yforehat = fn_forecast(Bhat,phil,nn);    % nfqm-by-nvar, in log
    else
        Xfexoe = fn_dataext(fdates(1,:),fdates(numel(fdates),:),xdatae(:,[1:2 2+nvar+1:2+nvar+nexo-1]));
        %Xfexoe = fn_dataext(fdates(1,:),fdates(end,:),xdatae(:,[1:2 2+nvar+1:2+nvar+nexo-1]));
        yforehat = fn_forecast(Bhat,phil,nn,nexo,Xfexoe(:,3:end));    % nfqm-by-nvar, in log
    end
    yforehate = [fdates yforehat];
    %
    yact1e = fn_dataext([yrEnd-nayr 1],[yrEnd qmEnd],xdatae(:,1:nvar+2));
    if options_.ms.real_pseudo_forecast
        %yact2e = fn_dataext([yrEnd-nayr 1],E2yrqm,xdatae);
        yact2e = fn_dataext([yrEnd-nayr 1],[fdates(end,1) options_.ms.freq],xdatae(:,1:nvar+2));
    else
        yact2e=yact1e;
    end
    yafhate = [yact1e; yforehate];  % actual and forecast
                                    %
                                    %===== Converted to mg, qg, and calendar yg
                                    %
    [yafyrghate,yafyrhate,yafqmyghate] = fn_datana(yafhate,options_.ms.freq,options_.ms.log_var(1:nlogeno),options_.ms.percent_var(1:npereno));
    % actual and forecast growth rates
    [yact2yrge,yact2yre,yact2qmyge] = fn_datana(yact2e,options_.ms.freq,options_.ms.log_var(1:nlogeno),options_.ms.percent_var(1:npereno));
    % only actual growth rates
    yafyrghate
    if options_.ms.indxgforhat
        keyindx = [1:nvar];
        conlab=['unconditional'];

        figure
        yafyrghate(:,3:end) = yafyrghate(:,3:end)/100;
        yact2yrge(:,3:end) = yact2yrge(:,3:end)/100;
        fn_foregraph(yafyrghate,yact2yrge,keyindx,rnum,cnum,options_.ms.freq,ylab,forelabel,conlab)
    end

    %-------------------------------------------------
    % Setup for point conditional forecast
    % ML Conditional Forecast
    %-------------------------------------------------
    %
    %% See Zha's note "Forecast (1)" p. 5, RATS manual (some errors in RATS), etc.
    %
    %% Some notations:  y(t+1) = y(t)B1 + e(t+1)inv(A0). e(t+1) is 1-by-n.
    %%    Let r(t+1)=e(t+1)inv(A0) + e(t+2)C + .... where inv(A0) is impulse
    %%          response at t=1, C at t=2, etc. The row of inv(A0) or C is
    %%          all responses to one shock.
    %%    Let r be q-by-1 (such as r(1) = r(t+1)
    %%                 = y(t+1) (constrained) - y(t+1) (forecast)).
    %%    Use impulse responses to find out R (k-by-q) where k=nvar*nsteps
    %%        where nsteps the largest constrained step.  The key of the program
    %%        is to creat R using impulse responses
    %%    Optimal solution for shock e where R'*e=r and e is k-by-1 is
    %%                 e = R*inv(R'*R)*r.
    %

    if (nconstr > 0)
        %*** initializing
        stepcon=cell(nconstr,1);  % initializing, value y conditioned
        valuecon=zeros(nconstr,1);  % initializing, value y conditioned
        varcon=zeros(nconstr,1);  % initializing, endogous variables conditioned
        varcon(:)=PorR;     % 1: Pcm; 3: FFR; 5: CPI

        %
        for i=1:nconstr
            if i<=nconstr1
                stepcon{i}=i;      % FFR
            elseif i<=2*nconstr1
                stepcon{i}=i-nconstr1;      % FFR
            elseif i<=3*nconstr1
                stepcon{i}=i-2*nconstr1;      % FFR
            elseif i<=4*nconstr1
                stepcon{i}=i-3*nconstr1;      % FFR
            elseif i<=5*nconstr1
                stepcon{i}=i-4*nconstr1;      % FFR
            elseif i<=6*nconstr1
                stepcon{i}=i-5*nconstr1;      % FFR
            end
        end

        %      for i=1:nconstr
        %         stepcon{i}=i;      % FFR
        %      end

        %      bend=12;
        %      stepcon{1}=[1:bend]'; % average over
        %      stepcon{nconstr1+1}=[1:options_.ms.freq-qmSub]';  % average over the remaing months in 1st forecast year
        %      stepcon{nconstr1+2}=[options_.ms.freq-qmSub+1:options_.ms.freq-qmSub+12]';    % average over 12 months next year
        %      stepcon{nconstr1+3}=[options_.ms.freq-qmSub+13:options_.ms.freq-qmSub+24]';    % average over 12 months. 3rd year
        %      stepcon{nconstr1+4}=[options_.ms.freq-qmSub+25:options_.ms.freq-qmSub+36]';    % average over 12 months. 4th year

        %      %**** avearage condition over, say, options_.ms.freq periods
        %      if qmEnd==options_.ms.freq
        %         stepcon{1}=[1:options_.ms.freq]';  % average over the remaing periods in 1st forecast year
        %      else
        %         stepcon{1}=[1:options_.ms.freq-qmEnd]';  % average over the remaing periods in 1st forecast year
        %      end
        %      for kj=2:nconstr
        %         stepcon{kj}=[length(stepcon{kj-1})+1:length(stepcon{kj-1})+options_.ms.freq]';    % average over 12 months next year
        %      end

        if options_.ms.real_pseudo_forecast
            %         %*** conditions in every period
            %         for i=1:nconstr
            %            valuecon(i) = yact(actup+i,varcon(i));
            %            %valuecon(i) = mean( yact(actup+1:actup+bend,varcon(i)) );
            %            %valuecon(i) = 0.060;      % 95:01
            %            %valuecon(i) = (0.0475+0.055)/2;   % 94:10
            %         end

            %         %*** average condtions over,say, options_.ms.freq periods.
            %         for i=nconstr1+1:nconstr1+nconstr2
            %            i=1;
            %            valuecon(nconstr1+i) = ( ( mean(ylast12Cal(:,varcon(nconstr1+i)),1) + ...
            %                 log(1+yactCalyg(yAg-yFg+i,varcon(nconstr1+i))/100) )*options_.ms.freq - ...
            %                 yCal_1(:,varcon(nconstr1+i)) ) ./ length(stepcon{nconstr1+i});
            %                             % the same as unconditional "yactCalyg" 1st calendar year
            %            i=2;
            %            valuecon(nconstr1+i) = mean(ylast12Cal(:,varcon(nconstr1+i))) +  ...
            %                 log(1+yactCalyg(yAg-yFg+1,varcon(nconstr1+i))/100) ...
            %                                + log(1+yactCalyg(yAg-yFg+i,varcon(nconstr1+i))/100);
            %                                    % the same as actual "yactCalgy" 2nd calendar year
            %            i=3;
            %            valuecon(nconstr1+i) = valuecon(nconstr1+i-1) + ...
            %                                        log(1+yactCalyg(yAg-yFg+i,varcon(nconstr1+i))/100);
            %                                    % the same as actual "yactCalgy" 3rd calendar year
            %            %i=4;
            %            %valuecon(nconstr1+i) = valuecon(nconstr1+i-1) + ...
            %            %                            log(1+yactCalyg(yAg-yFg+i,varcon(nconstr1+i))/100);
            %                                    % the same as actual "yactCalgy" 4th calendar year
            %         end

            %*** conditions in every period
            vpntM = fn_dataext(E1yrqm, E2yrqm,xdatae); % point value matrix with dates
                                                       %     vaveM = fn_dataext([yrEnd+1 0],[yrEnd+options_.forecast  0],yact2yre); % average value matrix with dates
            for i=1:nconstr
                if i<=nconstr1
                    valuecon(i) = vpntM(i,2+varcon(i)); % 2: first 2 elements are dates
                elseif i<=2*nconstr1
                    valuecon(i) = vpntM(i-nconstr1,2+varcon(i));
                elseif i<=3*nconstr1
                    valuecon(i) = vpntM(i-2*nconstr1,2+varcon(i));
                elseif i<=4*nconstr1
                    valuecon(i) = vpntM(i-3*nconstr1,2+varcon(i));
                elseif i<=5*nconstr1
                    valuecon(i) = vpntM(i-4*nconstr1,2+varcon(i));
                elseif i<=6*nconstr1
                    valuecon(i) = vpntM(i-5*nconstr1,2+varcon(i));
                end
            end

            %         %*** average condtions over,say, options_.ms.freq periods.
            %         if qmEnd==options_.ms.freq
            %            vaveM = fn_dataext([yrEnd+1 0],[yrEnd+options_.forecast  0],yact2yre); % average value matrix with dates
            %            valuecon(1) = vaveM(1,2+varcon(1));  % 2: first 2 elements are dates
            %         else
            %            vaveM = fn_dataext([yrEnd 0],[yrEnd+options_.forecast  0],yact2yre); % average value matrix with dates
            %            yactrem = fn_dataext([yrEnd qmEnd+1],[yrEnd options_.ms.freq],xdatae);
            %            valuecon(1) = sum(yactrem(:,2+varcon(1)),1)/length(stepcon{1});
            %                                    % 2: first 2 elements are dates
            %         end
            %         for kj=2:nconstr
            %            valuecon(kj) = vaveM(kj,2+varcon(kj));  % 2: first 2 elements are dates
            %         end
        else
            vpntM = dataext([yrEnd qmEnd+1],[yrEnd qmEnd+2],xdatae); % point value matrix with dates
            for i=1:nconstr
                if i<=nconstr1
                    valuecon(i) = vpntM(i,2+varcon(i)); % 2: first 2 elements are dates; Poil
                elseif i<=2*nconstr1
                    valuecon(i) = vpntM(i-nconstr1,2+varcon(i)); % 2: first 2 elements are dates; M2
                elseif i<=3*nconstr1
                    valuecon(i) = vpntM(i-2*nconstr1,2+varcon(i)); % 2: first 2 elements are dates; FFR
                elseif i<=4*nconstr1
                    valuecon(i) = vpntM(i-3*nconstr1,2+varcon(i)); % 2: first 2 elements are dates; CPI
                elseif i<=5*nconstr1
                    valuecon(i) = vpntM(i-4*nconstr1,2+varcon(i)); % 2: first 2 elements are dates; U
                elseif i<=5*nconstr1+nconstr2
                    valuecon(i)=xdata(end,5)+(i-5*nconstr1)*log(1.001)/options_.ms.freq;  %CPI
                elseif i<=5*nconstr1+2*nconstr2
                    valuecon(i)=0.0725;  %FFR
                else
                    valuecon(i)=xdata(end,6)+(i-5*nconstr1-2*nconstr2)*0.01/nfqm;  %U
                end
            end
            %valuecon(i) = 0.060;      % 95:01
        end
    else
        valuecon = [];
        stepcon = [];
        varcon = [];
    end

    nstepsm = 0;   % initializing, the maximum step in all constraints
    for i=1:nconstr
        nstepsm = max([nstepsm max(stepcon{i})]);
    end

    if cnum
        if options_.ms.real_pseudo_forecast && options_.ms.banact
            for i=1:length(banstp{1})
                banval{1}(1:length(banstp{1}),1) = ...
                    yactCalyg(yAg-yFg+1:yAg-yFg+length(banstp{1}),banvar(1)) - 2;
                banval{1}(1:length(banstp{1}),2) = ...
                    yactCalyg(yAg-yFg+1:yAg-yFg+length(banstp{1}),banvar(1)) + 2;
            end
        end
    end


    %===================================================
    % ML conditional forecast
    %===================================================
    %/*
    [ychat,Estr,rcon] = fn_fcstidcnd(valuecon,stepcon,varcon,nstepsm,...
                                     nconstr,options_.ms.eq_ms,nvar,options_.ms.nlags ,phil,0,0,yforehat,imf3shat,A0hat,Bhat,...
                                     nfqm,options_.ms.tlindx,options_.ms.tlnumber,options_.ms.ncms,options_.ms.eq_cms);
    ychate = [fdates ychat];
    yachate = [yact1e; ychate];  % actual and condtional forecast
                                 %===== Converted to mg, qg, and calendar yg
    [yacyrghate,yacyrhate,yacqmyghate] = fn_datana(yachate,options_.ms.freq,options_.ms.log_var(1:nlogeno),options_.ms.percent_var(1:npereno));
    % actual and conditional forecast growth rates
    if options_.ms.indxgdls && nconstr
        keyindx = [1:nvar];
        %  conlab=['conditional on' ylab{PorR(1)}];
        conlab=['v-conditions'];

        figure
        fn_foregraph(yafyrghate,yact2yrge,keyindx,rnum,cnum,options_.ms.freq,ylab,forelabel,conlab)
    end

    if options_.ms.ncsk
        Estr = zeros(nfqm,nvar);
        Estr(1:2,:) = [
            -2.1838      -1.5779      0.53064    -0.099425     -0.69269      -1.0391
            1.9407       3.3138     -0.10563     -0.55457     -0.68772       1.3534
                      ];
        Estr(3:6,3) = [0.5*ones(1,4)]';     % MD shocks

        Estr(3:10,2) = [1.5 1.5 1.5*ones(1,6)]';    % MS shocks

        %Estr(3:6,6) = 1*ones(4,1);      % U shocks
        %Estr(8:11,4) = 1*ones(4,1);      % y shocks

        %Estr(3:10,2) = [2.5 2.5 1.5*ones(1,6)]';    % MS shocks alone

        nn = [nvar options_.ms.noptions_.ms.nlags  nfqm];
        ycEhat = forefixe(A0hat,Bhat,phil,nn,Estr);
        ycEhate = [fdates ycEhat];
        yacEhate = [yact1e; ycEhate];  % actual and condtional forecast
                                       %===== Converted to mg, qg, and calendar yg
        [yacEyrghate,yacEyrhate,yacEqmyghate] = datana(yacEhate,options_.ms.freq,options_.ms.log_var(1:nlogeno),options_.ms.percent_var(1:npereno));
        % actual and conditional forecast growth rates
        disp([sprintf('%4.0f %2.0f %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f\n',yacEyrghate')])

        if 1
            keyindx = [1:nvar];
            %  conlab=['conditional on' ylab{PorR(1)}];
            conlab=['shock-conditions'];

            figure
            gyrfore(yacEyrghate,yact2yrge,keyindx,rnum,cnum,ylab,forelabel,conlab)
        end
    end

    %-----------------------------------------------------------
    % Compute structural shocks for the whole sample period excluding dummy observations.
    %-----------------------------------------------------------
    ywod = y(options_.ms.dummy_obs +1:end,:);     % without dummy observations
    phiwod=phi(options_.ms.dummy_obs +1:end,:);    % without dummy observations
    eplhat=ywod*A0hat-phiwod*Fhat;
    qmStartWod = mod(qmStart+options_.ms.nlags ,options_.ms.freq);
    if (~qmStartWod)
        qmStartWod = options_.ms.freq;
    end
    yrStartWod = yrStart + floor((qmStart+options_.ms.nlags -1)/options_.ms.freq);
    dateswod = fn_dataext([yrStartWod qmStartWod],[yrEnd qmEnd],xdatae(:,[1:2]));
    eplhate = [dateswod eplhat];

    Aphat = Fhat;
end

%----------------------------------------
% Tests for LR, HQ, Akaike, Schwarz.  The following gives a guidance.
%   But the computation has to be done in a different M file by exporting fhat's
%   from different idfile's.
%----------------------------------------
%
%if ~options_.ms.contemp_reduced_form
%   SpHR=A0in'*A0in;
%end
%%
%if ~isnan(SpHR) && ~options_.ms.contemp_reduced_form
%   warning(' ')
%   disp('Make sure you run the program with options_.ms.contemp_reduced_form =1 first.')
%   disp('Otherwise, the following test results such as Schwartz are incorrect.')
%   disp('All other results such as A0ml and imfs, however, are correct.')
%   disp('Press anykey to contintue or ctrl-c to stop now')
%   pause

%   load SpHUout

%   logLHU=-fss*sum(log(diag(chol(SpHU)))) -0.5*fss*nvar       % unrestricted logLH

%   logLHR=-fhat                                % restricted logLH
%   tra = reshape(SpHU,nvar*nvar,1)'*reshape(A0*A0',nvar*nvar,1);
%   df=(nvar*(nvar+1)/2 - length(a0indx));
%   S=2*(logLHU-logLHR);
%   SC = (nvar*(nvar+1)/2 - length(a0indx)) * log(fss);
%   disp(['T -- effective sample size:   ' num2str(fss)])
%   disp(['Trace in the overidentified posterior:   ' num2str(tra)])
%   disp(['Chi2 term -- 2*(logLHU-logLHR):   ' num2str(S)])
%   disp(['Degrees of freedom:   ' num2str(df)])
%   disp(['SC -- df*log(T):   ' num2str(SC)])
%   disp(['Akaike -- 2*df:   ' num2str(2*df)])
%   disp(['Classical Asymptotic Prob at chi2 term:   ' num2str(cdf('chi2',S,df))])

%   %*** The following is the eigenanalysis in the difference between
%   %***    unrestricted (U) and restricted (R)
%   norm(A0'*SpHU*A0-diag(diag(ones(6))))/6;
%   norm(SpHU-A0in'*A0in)/6;

%   corU = corr(SpHU);
%   corR = corr(SpHR);

%   [vU,dU]=eigsort(SpHU,1);
%   [vR,dR]=eigsort(SpHR,1);

%   [log(diag(dU)) log(diag(dR)) log(diag(dU))-log(diag(dR))];

%   sum(log(diag(dU)));
%   sum(log(diag(dR)));
%else
%   disp('To run SC test, turn options_.ms.contemp_reduced_form =1 first and then turn options_.ms.contemp_reduced_form =0')
%end


%***** Simply regression
%X=[phi(:,3) y(:,2)-phi(:,2) y(:,1)-phi(:,7) ones(fss,1)];
%� Y=y(:,3);
%� b=regress(Y,X)

%=== Computes the roots for the whole system.
rootsinv = fn_varoots(Bhat,nvar,options_.ms.nlags )
abs(rootsinv)


bhat =xhat;
n0const=n0;  % For constant parameter models.
n0const=n0;  % For constant parameter models.
npconst=np;  % For constant parameter models.
save outdata_a0dp_const.mat A0hat bhat Aphat n0const
