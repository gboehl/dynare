function [yy, xdir, isig, lam]=log_transform(y0,xdir0,isig,lam)
% [yy, xdir, isig, lam]=log_transform(y0,xdir0,isig,lam)
% Conduct automatic log transformation lam(yy/isig+lam)
% Inputs:
%   - y0    [double]    series to transform
%   - xdir  [char]      string indating the type of transformation:
%                           - log: standard log transformation
%                           - minuslog: log of minus (y0)
%                           - logsquared: log of y0^2
%                           - logskew: log of y0 shifted by lam
%   - isig  [double]    scaling factor for y0
%   - lam   [double]    shifting for y0
%
% Outputs:
%   - yy    [double]    transformed series
%   - xdir  [char]      string indating the type of transformation:
%                           - log: standard log transformation
%                           - minuslog: log of minus (y0)
%                           - logsquared: log of y0^2
%                           - logskew: log of y0 shifted by lam
%   - isig  [double]    scaling factor for y0
%   - lam   [double]    shifting for y0
%
% Notes: takes either one or four arguments. For one argument, the log
% transformation is conducted. For four arguments, the inverse
% transformation is applied.

% Written by Marco Ratto
% Joint Research Centre, The European Commission,
% marco.ratto@ec.europa.eu

% Copyright © 2012 European Commission
% Copyright © 2012-2017 Dynare Team
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

if nargin==4
    % inverse transformation
    yy = (exp(y0)-lam)*isig;
    return
end

if nargin==1
    xdir0='';
end
f=@(lam,y)gsa.skewness(log(y+lam));
isig=1;
if ~(max(y0)<0 || min(y0)>0)
    if gsa.skewness(y0)<0
        isig=-1;
        y0=-y0;
    end
    if isoctave
        n=hist(y0,10);
    else
        n=histcounts(y0,10);
    end
    if n(1)>20*n(end)
        try
            lam=fzero(f,[-min(y0)+10*eps -min(y0)+abs(median(y0))],[],y0);
        catch
            yl(1)=f(-min(y0)+10*eps,y0);
            yl(2)=f(-min(y0)+abs(median(y0)),y0);
            if abs(yl(1))<abs(yl(2))
                lam=-min(y0)+eps;
            else
                lam = -min(y0)+abs(median(y0));
            end
        end
        yy = log(y0+lam);
        xdir=[xdir0,'_logskew'];
    else
        isig=0;
        lam=0;
        yy = log(y0.^2);
        xdir=[xdir0,'_logsquared'];
    end
else
    if max(y0)<0
        isig=-1;
        y0=-y0;
        xdir=[xdir0,'_minuslog'];
    elseif min(y0)>0
        xdir=[xdir0,'_log'];
    end
    try
        lam=fzero(f,[-min(y0)+10*eps -min(y0)+median(y0)],[],y0);
    catch
        yl(1)=f(-min(y0)+10*eps,y0);
        yl(2)=f(-min(y0)+abs(median(y0)),y0);
        if abs(yl(1))<abs(yl(2))
            lam=-min(y0)+eps;
        else
            lam = -min(y0)+abs(median(y0));
        end
    end
    lam = max(lam,0);
    yy = log(y0+lam);
end
