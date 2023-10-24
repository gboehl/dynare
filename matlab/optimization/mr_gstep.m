function [f0, x, ig] = mr_gstep(h1,x,bounds,func0,penalty,htol0,Verbose,Save_files,gradient_epsilon,parameter_names,varargin)
% [f0, x, ig] = mr_gstep(h1,x,bounds,func0,penalty,htol0,Verbose,Save_files,gradient_epsilon,parameter_names,varargin)
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

% Copyright © 2006-2020 Dynare Team
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
            [f0, x, ~, ~] = csminit1(func0,x,penalty,f0,gg,0,diag(hh),Verbose,varargin{:});
            ig(i)=1;
            if Verbose
                fprintf(['Done for param %s = %8.4f\n'],parameter_names{i},x(i))
            end
        end
        xh1=x;
    end
    x = check_bounds(x,bounds);
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
    x(inx) = bounds(inx,2)-eps;
end

inx = find(x<=bounds(:,1));
if ~isempty(inx)
    x(inx) = bounds(inx,1)+eps;
end
