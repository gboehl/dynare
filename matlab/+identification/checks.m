function [condX, rankX, ind0, indno, ixno, Mco, Pco, jweak, jweak_pair] = checks(X, test_flag, tol_rank, tol_sv, param_nbr)
% function [condX, rankX, ind0, indno, ixno, Mco, Pco, jweak, jweak_pair] = checks(X, test_flag, tol_rank, tol_sv, param_nbr)
% -------------------------------------------------------------------------
% Checks rank criteria of identification tests and finds out parameter sets
% that are not identifiable via the nullspace, pairwise correlation
% coefficients and multicorrelation coefficients
% =========================================================================
% INPUTS
%    * X               [matrix] dependent on test_flag:
%                               test_flag = 0: Sample information matrix (Ahess)
%                               test_flag = 1: Jacobian of Moments (J), reduced-form (dTAU) or dynamic model (dLRE)
%                               test_flag = 2: Jacobian of minimal system (D)
%                               test_flag = 3: Gram matrix (hessian or correlation type matrix) of spectrum (G)
% -------------------------------------------------------------------------
% OUTPUTS
%    * cond             [double]        condition number of X
%    * rank             [double]        rank of X with tolerance tol_rank
%    * ind0             [vector]        binary indicator for non-zero columns of H
%    * indno            [matrix]        index of non-identified params
%    * ixno             [integer]       number of rows in indno
%    * Mco              [matrix]        multicollinearity coefficients
%    * Pco              [matrix]        pairwise correlations
%    * jweak            [matrix]        gives 1 if the parameter has Mco=1 (with tolerance tol_rank)
%    * jweak_pair       [(vech) matrix] gives 1 if a couple parameters has Pco=1 (with tolerance tol_rank)
% -------------------------------------------------------------------------
% This function is called by
%   * identification.analysis.m
% -------------------------------------------------------------------------
% This function calls
%    * identification.cosn
%    * dyn_vech
%    * vnorm
% =========================================================================
% Copyright © 2010-2019 Dynare Team
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
% =========================================================================
if issparse(X)
    X = full(X);
end
if nargin < 3 || isempty(tol_rank) || strcmp(tol_rank,'robust')
    tol_rank = max(size(X)) * eps(norm(X)); %tolerance level used for rank computations
end
if nargin < 4 || isempty(tol_sv)
    tol_sv = 1.e-3; %tolerance level for zero singular value
end
if nargin < 5 || isempty(param_nbr)
    param_nbr = size(X,2);
end

indno = zeros(1,param_nbr);
rankrequired = size(X,2);

if test_flag == 2
    % Komunjer and Ng's D
    Xpar = X(:,1:param_nbr);
    Xrest = X(:,(param_nbr+1):end);
else
    Xpar = X;
    Xrest = [];
end

% find non-zero columns at machine precision
if size(Xpar,1) > 1
    ind1 = find(identification.vnorm(Xpar) >= eps);
else
    ind1 = find(abs(Xpar) >= eps);   % if only one parameter
end

if test_flag == 3
    Xparnonzero = Xpar(ind1,ind1); % focus on non-zero rows and columns for Qu and Tkachenko's G
else
    Xparnonzero = Xpar(:,ind1); % focus on non-zero columns
end

[~, ~, ee1] = svd( [Xparnonzero Xrest], 0 );
condX = cond([Xparnonzero Xrest]);
rankX = rank(X, tol_rank);
icheck = 0; %initialize flag which is equal to 0 if we already found all single parameters that are not identified
if param_nbr > 0 && (rankX<rankrequired)
    % search for singular values associated to ONE individual parameter
    % Compute an orthonormal basis for the null space using the columns of ee1 that correspond
    % to singular values equal to zero and associated to an individual parameter
    ee0 = rankX+1:size([Xparnonzero Xrest],2); %look into last columns with singular values of problematic parameter sets (except single parameters)
    ind11 = ones(length(ind1),1); %initialize
    for j=1:length(ee0)
        % check if nullspace is spanned by only one parameter
        if length(find(abs(ee1(:,ee0(j))) > tol_sv))==1 %note that tol_sv is not the same tolerance used for rank computations
            icheck=1; %indicate that additional single parameters are found
            if test_flag == 2
                temp = (abs(ee1(:,ee0(j))) <= tol_sv);
                ind11 = ind11.*temp(1:(end-size(Xrest,2))); % find non-zero columns
            else
                ind11 = ind11.*(abs(ee1(:,ee0(j))) <= tol_sv); % find non-zero columns
            end
        end
    end
    ind1 = ind1(find(ind11)); % find non-zero columns
end

if icheck
    %if we found additional single parameters we need to recheck if we found all parameters
    if test_flag == 3
        Xparnonzero = Xpar(ind1,ind1); % focus on non-zero rows and columns for Qu and Tkachenko's G
    else
        Xparnonzero = Xpar(:,ind1); % focus on non-zero columns
    end
    [~, ~, ee1] = svd( [Xparnonzero Xrest], 0 );
    condX = cond([Xparnonzero Xrest]);
    rankX = rank(X,tol_rank);
end

ind0 = zeros(1,param_nbr); %initialize
ind0(ind1) = 1;

% find near linear dependence problems via multicorrelation coefficients
if test_flag == 0 || test_flag == 3 % G is a Gram matrix and hence should be a correlation-like matrix
    if test_flag == 3 % For Qu and Tkachenko's G matrix we need to keep track of all parameters
        Mco = NaN(param_nbr,1);
    end
    Pco=NaN(param_nbr,param_nbr);         % pairwise correlation coefficient
    deltaX = sqrt(diag(X(ind1,ind1)));
    tildaX = X(ind1,ind1)./((deltaX)*(deltaX'));
    Mco(ind1,1)=(1-1./diag(inv(tildaX))); % multicorrelation coefficent
    Pco(ind1,ind1)=inv(X(ind1,ind1));
    sd=sqrt(diag(Pco));
    Pco = abs(Pco./((sd)*(sd')));
else
    Mco = NaN(param_nbr,1);
    for ii = 1:size(Xparnonzero,2)
        Mco(ind1(ii),:) = identification.cosn([Xparnonzero(:,ii) , Xparnonzero(:,find([1:1:size(Xparnonzero,2)]~=ii)), Xrest]);
    end
end

%% find out which parameters are involved
ixno = 0; %initialize number of rows
if param_nbr>0 && (rankX<rankrequired || min(1-Mco)<tol_rank)
    if length(ind1)<param_nbr
        % single parameters with zero columns
        ixno = ixno + 1;
        indno(ixno,:) = (~ismember(1:param_nbr,ind1));
    end
    ee0 = rankX+1:size([Xparnonzero Xrest],2); %look into last columns with singular values of problematic parameter sets (except single parameters)
    for j=1:length(ee0)
        % linearly dependent parameters
        ixno = ixno + 1;
        if test_flag == 2
            temp = (abs(ee1(:,ee0(j))) > tol_sv)';
            indno(ixno,ind1) = temp(1:(end-size(Xrest,2)));
        else
            indno(ixno,ind1) = (abs(ee1(:,ee0(j))) > tol_sv)';
        end
    end
end

%% here there is no exact linear dependence, but there are several near-dependencies, mostly due to strong pairwise colliniearities
jweak      = zeros(1,param_nbr);
jweak_pair = zeros(param_nbr,param_nbr);

if test_flag ~= 0
    % these tests only apply to Jacobians, not to Gram matrices, i.e. Hessian-type or 'covariance' matrices
    Pco = NaN(param_nbr,param_nbr);
    for ii = 1:size(Xparnonzero,2)
        Pco(ind1(ii),ind1(ii)) = 1;
        for jj = ii+1:size(Xparnonzero,2)
            Pco(ind1(ii),ind1(jj)) = identification.cosn([Xparnonzero(:,ii),Xparnonzero(:,jj),Xrest]);
            Pco(ind1(jj),ind1(ii)) = Pco(ind1(ii),ind1(jj));
        end
    end

    for j=1:param_nbr
        if Mco(j)>(1-tol_rank)
            jweak(j)=1;
            [~, jpair] = find(Pco(j,j+1:end)>(1-tol_rank));
            for jx=1:length(jpair)
                jweak_pair(j, jpair(jx)+j)=1;
                jweak_pair(jpair(jx)+j, j)=1;
            end
        end
    end
end

jweak_pair=dyn_vech(jweak_pair)'; % focus only on unique combinations
