function [hessian_mat, gg, htol1, ihh, hh_mat0, hh1, hess_info] = mr_hessian(x,func,penalty,hflag,htol0,hess_info,bounds,prior_std,Save_files,varargin)
% function [hessian_mat, gg, htol1, ihh, hh_mat0, hh1, hess_info] = mr_hessian(x,func,penalty,hflag,htol0,hess_info,bounds,prior_std,Save_files,varargin)
%  numerical gradient and Hessian, with 'automatic' check of numerical
%  error
%
% adapted from Michel Juillard original routine hessian.m
%
% Inputs:
%  - x                  parameter values
%  - func               function handle. The function must give two outputs:
%                       the log-likelihood AND the single contributions at times t=1,...,T
%                       of the log-likelihood to compute outer product gradient
%  - penalty            penalty due to error code
%  - hflag              0: Hessian computed with outer product gradient, one point
%                           increments for partial derivatives in gradients
%                       1: 'mixed' Hessian: diagonal elements computed with numerical second order derivatives
%                           with correlation structure as from outer product gradient;
%                           two point evaluation of derivatives for partial derivatives
%                           in gradients
%                       2: full numerical Hessian, computes second order partial derivatives
%                           uses Abramowitz and Stegun (1965) formulas 25.3.24 and 25.3.27
%                           p. 884.
%  - htol0              'precision' of increment of function values for numerical
%                       derivatives
%  - hess_info              structure storing the step sizes for
%                           computation of Hessian
%  - bounds                 prior bounds of parameters
%  - prior_std              prior standard devation of parameters (can be NaN)
%  - Save_files             indicator whether files should be saved
%  - varargin               other inputs
%                           e.g. in dsge_likelihood
%                           varargin{1} --> dataset_
%                           varargin{2} --> dataset_info
%                           varargin{3} --> options_
%                           varargin{4} --> M_
%                           varargin{5} --> estim_params_
%                           varargin{6} --> bayestopt_
%                           varargin{7} --> BoundsInfo
%                           varargin{8} --> oo_
%
% Outputs
%  - hessian_mat        hessian
%  - gg                 Jacobian
%  - htol1              updated 'precision' of increment of function values for numerical
%                       derivatives
%  - ihh                inverse outer product with modified std's
%  - hh_mat0            outer product hessian with modified std's
%  - hh1                updated hess_info.h1
%  - hess_info          structure with updated step length

% Copyright © 2004-2017 Dynare Team
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

n=size(x,1);

[f0,~, ff0]=penalty_objective_function(x,func,penalty,varargin{:});
h2=bounds(:,2)-bounds(:,1);
hmax=bounds(:,2)-x;
hmax=min(hmax,x-bounds(:,1));
if isempty(ff0)
    outer_product_gradient=0;
else
    outer_product_gradient=1;
end


hess_info.h1 = min(hess_info.h1,0.9.*hmax);

if htol0<hess_info.htol
    hess_info.htol=htol0;
end
if not(isnan(hess_info.htol_rel))
    htol1=hess_info.htol;
    hess_info.htol=abs(hess_info.htol_rel*f0);
end

xh1=x;
f1=zeros(size(f0,1),n);
f_1=f1;
if outer_product_gradient
    ff1=zeros(size(ff0));
    ff_1=ff1;
    ggh=zeros(size(ff0,1),n);
end

i=0;
hhtol=hess_info.htol*ones(n,1);
while i<n
    i=i+1;
    hess_info.htol=hhtol(i);
    h10=hess_info.h1(i);
    hcheck=0;
    xh1(i)=x(i)+hess_info.h1(i);
    [fx,~,ffx]=penalty_objective_function(xh1,func,penalty,varargin{:});
    it=1;
    dx=(fx-f0);
    ic=0;
    icount = 0;
    h0=hess_info.h1(i);
    istoobig=false;
    istoosmall=false;
    give_it_up=false;
    while (abs(dx(it))<0.5*hess_info.htol || abs(dx(it))>(3*hess_info.htol)) && icount<10 && ic==0 && hmax(i)>2*1.e-10 && not(give_it_up)
        icount=icount+1;
        istoobig(it)=false;
        istoosmall(it)=false;
        if abs(dx(it))<0.5*hess_info.htol
            istoosmall(it)=true;
            if hess_info.h1(i)==0.9*hmax(i)% || hess_info.h1(i)==0.3*abs(x(i))
                give_it_up=true;
            end
            if abs(dx(it)) ~= 0
                %                 htmp=min(max(1.e-10,0.3*abs(x(i))), 0.9*hess_info.htol/abs(dx(it))*hess_info.h1(i));
                htmp=0.9*hess_info.htol/abs(dx(it))*hess_info.h1(i);
            else
                htmp=2.1*hess_info.h1(i);
            end
            htmp = min(htmp,0.9*hmax(i));
            if any(h0(istoobig(1:it))) %&& htmp>=min(h0(istoobig(1:it)))
                htmp = 0.5*(min(h0(istoobig))+max(h0(istoosmall)));
            end
            hess_info.h1(i) = max(htmp,1.e-10);
            xh1(i)=x(i)+hess_info.h1(i);
            [fx,~,ffx]=penalty_objective_function(xh1,func,penalty,varargin{:});
        end
        if abs(dx(it))>(3*hess_info.htol)
            istoobig(it)=true;
            htmp= hess_info.htol/abs(dx(it))*hess_info.h1(i);
            if any(h0(istoosmall(1:it))) %&& htmp<=max(h0(istoosmall(1:it)))
                htmp = 0.5*(min(h0(istoobig))+max(h0(istoosmall)));
            end
            hess_info.h1(i) = max(htmp,1e-10);
            xh1(i)=x(i)+hess_info.h1(i);
            [fx,~,ffx]=penalty_objective_function(xh1,func,penalty,varargin{:});
            iter=0;
            while (fx-f0)==0 && iter<50
                hess_info.h1(i)= hess_info.h1(i)*2;
                xh1(i)=x(i)+hess_info.h1(i);
                [fx,~,ffx]=penalty_objective_function(xh1,func,penalty,varargin{:});
                ic=1;
                iter=iter+1;
            end
        end
        it=it+1;
        dx(it)=(fx-f0);
        h0(it)=hess_info.h1(i);
        if (hess_info.h1(i)<1.e-12*min(1,h2(i)) && hess_info.h1(i)<0.9*hmax(i))
            ic=1;
            hcheck=1;
        end
    end
    if icount == 10 || hess_info.h1(i)==1.e-10
        istoobig(it)=false;
        istoosmall(it)=false;
        if abs(dx(it))<0.5*hess_info.htol
            istoosmall(it)=true;
        end
        if abs(dx(it))>(3*hess_info.htol)
            istoobig(it)=true;
        end
        if any(istoobig) && (istoosmall(it) || istoobig(it))
            % always better to be wrong from above, to avoid numerical noise
            [ddx, ij]=min(dx(istoobig));
            fx=f0+ddx;
            hess_info.h1(i)=h0(ij);
        end
    end
    f1(:,i)=fx;
    if outer_product_gradient
        if any(isnan(ffx)) || isempty(ffx)
            ff1=ones(size(ff0)).*fx/length(ff0);
        else
            ff1=ffx;
        end
    end
    xh1(i)=x(i)-hess_info.h1(i);
    [fx,~,ffx]=penalty_objective_function(xh1,func,penalty,varargin{:});
    f_1(:,i)=fx;
    if outer_product_gradient
        if any(isnan(ffx)) || isempty(ffx)
            ff_1=ones(size(ff0)).*fx/length(ff0);
        else
            ff_1=ffx;
        end
        ggh(:,i)=(ff1-ff_1)./(2.*hess_info.h1(i));
    end
    xh1(i)=x(i);
    if hcheck && hess_info.htol<1
        hess_info.htol=min(1,max(min(abs(dx))*2,hess_info.htol*10));
        hess_info.h1(i)=h10;
        hhtol(i) = hess_info.htol;
        i=i-1;
    end
end

h_1=hess_info.h1;
xh1=x;
xh_1=xh1;

gg=(f1'-f_1')./(2.*hess_info.h1);

if outer_product_gradient
    if hflag==2
        % full numerical Hessian
        hessian_mat = zeros(size(f0,1),n*n);
        for i=1:n
            if i > 1
                k=i:n:n*(i-1);
                hessian_mat(:,(i-1)*n+1:(i-1)*n+i-1)=hessian_mat(:,k);
            end
            hessian_mat(:,(i-1)*n+i)=(f1(:,i)+f_1(:,i)-2*f0)./(hess_info.h1(i)*h_1(i));
            temp=f1+f_1-f0*ones(1,n);
            for j=i+1:n
                xh1(i)=x(i)+hess_info.h1(i);
                xh1(j)=x(j)+h_1(j);
                xh_1(i)=x(i)-hess_info.h1(i);
                xh_1(j)=x(j)-h_1(j);
                temp1 = penalty_objective_function(xh1,func,penalty,varargin{:});
                temp2 = penalty_objective_function(xh_1,func,penalty,varargin{:});
                hessian_mat(:,(i-1)*n+j)=-(-temp1 -temp2+temp(:,i)+temp(:,j))./(2*hess_info.h1(i)*h_1(j));
                xh1(i)=x(i);
                xh1(j)=x(j);
                xh_1(i)=x(i);
                xh_1(j)=x(j);
            end
        end
    elseif hflag==1
        % full numerical 2nd order derivs only in diagonal
        hessian_mat = zeros(size(f0,1),n*n);
        for i=1:n
            dum = (f1(:,i)+f_1(:,i)-2*f0)./(hess_info.h1(i)*h_1(i));
            hessian_mat(:,(i-1)*n+i)=dum;
            if any(dum<=eps)
                hessian_mat(dum<=eps,(i-1)*n+i)=max(eps, gg(i)^2);
            end
        end
    end

    gga=ggh.*kron(ones(size(ff1)),2.*hess_info.h1');  % re-scaled gradient
    hh_mat=gga'*gga;  % rescaled outer product hessian
    hh_mat0=ggh'*ggh;  % outer product hessian
    A=diag(2.*hess_info.h1);  % rescaling matrix
    % igg=inv(hh_mat);  % inverted rescaled outer product hessian
    ihh=A'*(hh_mat\A);  % inverted outer product hessian (based on rescaling)
    if hflag>0 && min(eig(reshape(hessian_mat,n,n)))>0
        hh0 = A*reshape(hessian_mat,n,n)*A';  %rescaled second order derivatives
        hh = reshape(hessian_mat,n,n);  %second order derivatives
        sd0=sqrt(diag(hh0));   %rescaled 'standard errors' using second order derivatives
        sd=sqrt(diag(hh_mat));  %rescaled 'standard errors' using outer product
        hh_mat=hh_mat./(sd*sd').*(sd0*sd0');  %rescaled inverse outer product with 'true' std's
        ihh=A'*(hh_mat\A);  % update inverted outer product hessian with 'true' std's
        sd=sqrt(diag(ihh));   %standard errors
        sdh=sqrt(1./diag(hh));   %diagonal standard errors
        for j=1:length(sd)
            % some heuristic normalizations of the standard errors that
            % avoid numerical issues in outer product
            sd0(j,1)=min(prior_std(j), sd(j));  %prior std
            sd0(j,1)=10^(0.5*(log10(sd0(j,1))+log10(sdh(j,1))));
        end
        inv_A=inv(A);
        ihh=ihh./(sd*sd').*(sd0*sd0');  %inverse outer product with modified std's
        igg=inv_A'*ihh*inv_A;  % inverted rescaled outer product hessian with modified std's
        % hh_mat=inv(igg);   % outer product rescaled hessian with modified std's
        hh_mat0=inv_A'/igg*inv_A;  % outer product hessian with modified std's
        %     sd0=sqrt(1./diag(hh0));   %rescaled 'standard errors' using second order derivatives
        %     sd=sqrt(diag(igg));  %rescaled 'standard errors' using outer product
        %     igg=igg./(sd*sd').*(sd0*sd0');  %rescaled inverse outer product with 'true' std's
        %     hh_mat=inv(igg);   % rescaled outer product hessian with 'true' std's
        %     ihh=A'*igg*A;  % inverted outer product hessian
        %     hh_mat0=inv(A)'*hh_mat*inv(A);  % outer product hessian with 'true' std's
    end
    if hflag<2
        hessian_mat=hh_mat0(:);
    end

    if any(isnan(hessian_mat))
        hh_mat0=eye(length(hh_mat0));
        ihh=hh_mat0;
        hessian_mat=hh_mat0(:);
    end
    hh1=hess_info.h1;
    if Save_files
        save('hess.mat','hessian_mat')
    end
else
    hessian_mat=[];
    ihh=[];
    hh_mat0 = [];
    hh1 = [];
end

if isnan(hess_info.htol_rel)
    htol1=hhtol;
end
