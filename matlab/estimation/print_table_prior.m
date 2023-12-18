function print_table_prior(lb, ub, options_, M_, bayestopt_, estim_params_)
% print_table_prior(lb, ub, options_, M_, bayestopt_, estim_params_)
% This routine prints in the command window some descriptive statistics about the prior distribution.
% Inputs:
%  o lb             [double]    lower bound    
%  o ub             [double]    upper bound
%  o M_             [structure] Definition of the model
%  o bayestopt_     [structure] describing the priors
%  o estim_params_  [structure] characterizing parameters to be estimated

% Copyright © 2015-2023 Dynare Team
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

PriorNames = 'Beta';
PriorNames = strvcat(PriorNames,'Gamma');
PriorNames = strvcat(PriorNames,'Gaussian');
PriorNames = strvcat(PriorNames,'Inverted Gamma');
PriorNames = strvcat(PriorNames,'Uniform');
PriorNames = strvcat(PriorNames,'Inverted Gamma -- 2');
PriorNames = strvcat(PriorNames,'Dirichlet');
PriorNames = strvcat(PriorNames,'Weibull');

n = size(bayestopt_.name,1); % Numbe rof estimated parameters.

l1 = printline(10, '-');
T1 = strvcat(l1, 'PARAMETER ');
T1 = strvcat(T1, l1);

l2 = printline(133, '-');
T2 = strvcat(l2, sprintf('Prior shape      \t Prior mean \t Prior mode \t Prior std. \t Prior lb \t Prior ub \t Prior HPD lb \t Prior HPH ub'));
T2 = strvcat(T2, l2);

prior_trunc_backup = options_.prior_trunc ;
options_.prior_trunc = (1-options_.prior_interval)/2 ;
PriorIntervals = prior_bounds(bayestopt_, options_.prior_trunc) ;
options_.prior_trunc = prior_trunc_backup ;

RESIZE = false;

for i=1:size(bayestopt_.name,1)
    Name = get_the_name(i,1,M_,estim_params_,options_.varobs);
    if length(Name)>size(T1,2)
        resize = true;
    else
        resize = false;
    end
    T1 = strvcat(T1, Name);
    if resize
        RESIZE = true;
        l1 = printline(length(Name));
        T1(1,:) = l1;
        T1(3,:) = l1;
    end
    PriorShape = PriorNames(bayestopt_.pshape(i),:);
    PriorMean = bayestopt_.p1(i);
    PriorMode = bayestopt_.p5(i);
    PriorStandardDeviation = bayestopt_.p2(i);
    switch bayestopt_.pshape(i)
      case { 1 , 5 }
        LowerBound = bayestopt_.p3(i);
        UpperBound = bayestopt_.p4(i);
        if ~isinf(lb(i))
            LowerBound=max(LowerBound,lb(i));
        end
        if ~isinf(ub(i))
            UpperBound=min(UpperBound,ub(i));
        end
      case { 2 , 4 , 6 , 8}
        LowerBound = bayestopt_.p3(i);
        if ~isinf(lb(i))
            LowerBound=max(LowerBound,lb(i));
        end
        if ~isinf(ub(i))
            UpperBound=ub(i);
        else
            UpperBound = Inf;
        end
      case 3
        if isinf(bayestopt_.p3(i)) && isinf(lb(i))
            LowerBound = -Inf;
        else
            LowerBound = bayestopt_.p3(i);
            if ~isinf(lb(i))
                LowerBound=max(LowerBound,lb(i));
            end
        end
        if isinf(bayestopt_.p4(i)) && isinf(ub(i))
            UpperBound = Inf;
        else
            UpperBound = bayestopt_.p4(i);
            if ~isinf(ub(i))
                UpperBound=min(UpperBound,ub(i));
            end
        end
      otherwise
        error('get_prior_info:: Dynare bug!')
    end
    format_string = build_format_string(PriorMode, PriorStandardDeviation,LowerBound,UpperBound);
    str = sprintf(format_string, ...
                  PriorShape, ...
                  PriorMean, ...
                  PriorMode, ...
                  PriorStandardDeviation, ...
                  LowerBound, ...
                  UpperBound, ...
                  PriorIntervals.lb(i), ...
                  PriorIntervals.ub(i) );
    T2 = strvcat(T2, str);
end

T1 = strvcat(T1, l1);
T2 = strvcat(T2, l2);

skipline(2)

if RESIZE
    l0 = printline(2);
    T0 = strvcat(l0,'  ');
    T0 = strvcat(T0, l0);
    T0 = strvcat(T0, repmat('  ', n, 1));
    T0 = strvcat(T0, l0);
    disp([T1, T0, T2])
else
    disp([T1, T2])
end

skipline(2)


function format_string = build_format_string(PriorMode,PriorStandardDeviation,LowerBound,UpperBound)
format_string = ['%s \t %6.4f \t'];
if isnan(PriorMode)
    format_string = [ format_string , ' %s \t'];
else
    format_string = [ format_string , ' %6.4f \t'];
end
if ~isnumeric(PriorStandardDeviation)
    format_string = [ format_string , ' %s \t'];
else
    format_string = [ format_string , ' %6.4f \t'];
end
if ~isnumeric(LowerBound)
    format_string = [ format_string , ' %s \t'];
else
    format_string = [ format_string , ' %6.4f \t'];
end
if ~isnumeric(UpperBound)
    format_string = [ format_string , ' %s \t'];
else
    format_string = [ format_string , ' %6.4f \t'];
end
format_string = [ format_string , ' %6.4f \t %6.4f'];
