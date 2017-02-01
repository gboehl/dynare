function [ya, yass, gya, gyass] = quarterly2annual(y,yss,GYTREND0,type,islog)
% function [ya, yass, gya, gyass] = quarterly2annual(y,yss,GYTREND0,type,islog)
% transforms quarterly (log-)level time series to annual level and growth rate
% it accounts for stock/flow/deflator series.
%
% INPUTS
% y        quarterly time series
% yss      steady state of y
% GYTREND0 growth rate of y
% type     1 sum (default)
%          2 average
%          3 last period (Q4)
%          4 geometric average
% islog    0 level (default)
%          1 log-level
%
% OUTPUTS
% ya       annual (log-)level
% yass     annual steadystate (log-)level  
% gya      annual growth rate
% gyass    annual growth rate steadystate

% Copyright (C) 2017 Dynare Team
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

if nargin ==0,
    disp('[ya, yass, gya, gyass] = quarterly2annual(y,yss,GYTREND0,type,islog);')
    return
end

if nargin<4 || isempty(type),
    type=1;
end
if nargin<5 || isempty(islog),
    islog=0;
end
if islog
    y=exp(y+yss);
    yss=exp(yss);
    y=y-yss;
end
switch type
    case 1
        yass = yss*(exp(-GYTREND0*3)+exp(-GYTREND0*2)+exp(-GYTREND0)+1);
        tmp = lagged(y,3)*exp(-GYTREND0*3)+lagged(y,2)*exp(-GYTREND0*2)+lagged(y,1)*exp(-GYTREND0)+y; % annualized level
    case 2
        yass = yss*(exp(-GYTREND0*3)+exp(-GYTREND0*2)+exp(-GYTREND0)+1)/4;
        tmp = (lagged(y,3)*exp(-GYTREND0*3)+lagged(y,2)*exp(-GYTREND0*2)+lagged(y,1)*exp(-GYTREND0)+y)/4; % annualized level
    case 3
        yass=yss;
        tmp = y;
    case 4
        yass = yss*(exp(-GYTREND0*3/2));
        tmp = (lagged(y+yss,3)*exp(-GYTREND0*3).*lagged(y+yss,2)*exp(-GYTREND0*2).*lagged(y+yss,1)*exp(-GYTREND0).*(y+yss)).^(1/4); % annualized level        
        tmp = tmp - yass;
    otherwise
        error('Wrong type input')
end

ya = tmp(4:4:end);
% annual growth rate
gyass = GYTREND0*4;
gya = (ya+yass)./(lagged(ya,1)+yass).*exp(4*GYTREND0)-1-gyass;

if islog
    ya=log(ya+yass);
    yass=log(yass);
    ya=ya-yass;
end

