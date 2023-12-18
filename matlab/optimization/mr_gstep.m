function [f0, x, ig] = mr_gstep(h1,x,bounds,func0,penalty,htol0,Verbose,Save_files,gradient_epsilon,parameter_names,robust,varargin)
% [f0, x, ig] = mr_gstep(h1,x,bounds,func0,penalty,htol0,Verbose,Save_files,gradient_epsilon,parameter_names,robust,varargin)
%
% Gibbs type step in optimisation
%
% varargin{1} --> dataset_
% varargin{2} --> dataset_info
% varargin{3} --> options_
% varargin{4} --> M_
% varargin{5} --> estim_params_
% varargin{6} --> bayestopt_
% varargin{7} --> BoundsInfo
% varargin{8} --> oo_

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

n=size(x,1);
if isempty(h1)
    h1=gradient_epsilon*ones(n,1);
end


if isempty(htol0)
    htol = 1.e-6;
else
    htol = htol0;
end
if length(htol)==1
    htol=htol*ones(n,1);
end
f0=penalty_objective_function(x,func0,penalty,varargin{:});

xh1=x;
f1=zeros(size(f0,1),n);
f_1=f1;

i=0;
ig=zeros(n,1);
while i<n
    i=i+1;
    h10=h1(i);
    hcheck=0;
    dx=[];
    xh1(i)=x(i)+h1(i);
    fx = penalty_objective_function(xh1,func0,penalty,varargin{:});
    f1(:,i)=fx;
    xh1(i)=x(i)-h1(i);
    fx = penalty_objective_function(xh1,func0,penalty,varargin{:});
    f_1(:,i)=fx;
    if hcheck && htol(i)<1
        htol(i)=min(1,max(min(abs(dx))*2,htol(i)*10));
        h1(i)=h10;
        xh1(i)=x(i);
        i=i-1;
    else
        gg=zeros(size(x));
        hh=gg;
        gg(i)=(f1(i)'-f_1(i)')./(2.*h1(i));
        hh(i) = 1/max(1.e-9,abs( (f1(i)+f_1(i)-2*f0)./(h1(i)*h1(i)) ));
        if gg(i)*(hh(i)*gg(i))/2 > htol(i)
            [ff, xx,~,retcode] = csminit1(func0,x,penalty,f0,gg,0,diag(hh),Verbose,varargin{:});
            if retcode && robust 
                if abs(x(i))<1.e-6
                    xa=transpose(linspace(x(i)/2, sign(x(i))*1.e-6*3/2, 7));
                else
                    xa=transpose(linspace(x(i)/2, x(i)*3/2, 7));
                end
                fa=NaN(7,1);
                for k=1:7
                    xh1(i)=xa(k);
                    fa(k,1) = penalty_objective_function(xh1,func0,penalty,varargin{:});
                end
                b=[ones(7,1) xa xa.*xa./2]\fa;
                gg(i)=x(i)*b(3)+b(2);
                hh(i)=1/b(3);
                [ff2, xx2] = csminit1(func0,x,penalty,f0,gg,0,diag(hh),Verbose,varargin{:});
                if ff2<ff
                    ff=ff2;
                    xx=xx2;
                end
                if min(fa)<ff
                    [ff, im]=min(fa);
                    xx(i)=xa(im);
                end
            end
            ig(i)=1;
            if robust
            if not(isequal(xx , check_bounds(xx,bounds)))
                xx = check_bounds(xx,bounds);
                if xx(i)<x(i)   
                    % lower bound
                    xx(i) = min(xx(i)+h1(i), 0.5*(xx(i)+x(i)));
                else
                    % upper bound
                    xx(i) = max(xx(i)-h1(i), 0.5*(xx(i)+x(i)));
                end
                [ff,exit_flag]=penalty_objective_function(xx,func0,penalty,varargin{:});
                if exit_flag~=1
                    disp_verbose('last step exited with bad status!',Verbose)
                elseif ff<f0
                    f0=ff;
                    x=xx;
                end
            else
                % check improvement wrt predicted one
                if abs(f0-ff) < abs(gg(i)*(hh(i)*gg(i))/2/100) || abs(x(i)-xx(i))<1.e-10
                    [ff1, xx1] = csminit1(func0,x,penalty,f0,-gg,0,diag(hh),Verbose,varargin{:});
                    if not(isequal(xx1 , check_bounds(xx1,bounds)))
                        xx1 = check_bounds(xx1,bounds);
                        if xx1(i)<x(i)
                            % lower bound
                            xx1(i) = min(xx1(i)+h1(i), 0.5*(xx1(i)+x(i)));
                        else
                            % upper bound
                            xx1(i) = max(xx1(i)-h1(i), 0.5*(xx1(i)+x(i)));
                        end
                        [ff1,exit_flag]=penalty_objective_function(xx1,func0,penalty,varargin{:});
                        if exit_flag~=1
                            disp_verbose('last step exited with bad status!',Verbose)
                        end
                    end
                    if ff1<ff
                        ff=ff1;
                        xx=xx1;
                    end
                end
                f0=ff;
                x=xx;
            end
            else
                f0=ff;
                x=xx;
                x = check_bounds(x,bounds);
            end
            if Verbose
                fprintf('Done for param %s = %8.4f; f = %8.4f\n',parameter_names{i},x(i),f0)
            end
        end
        xh1=x;
    end
    if Save_files
        save('gstep.mat','x','h1','f0')
    end
end
if Save_files
    save('gstep.mat','x','h1','f0')
end

return


function x = check_bounds(x,bounds)

inx = find(x>=bounds(:,2));
if ~isempty(inx)
    x(inx) = bounds(inx,2)-1.e-10;
end

inx = find(x<=bounds(:,1));
if ~isempty(inx)
    x(inx) = bounds(inx,1)+1.e-10;
end
