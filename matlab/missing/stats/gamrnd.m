function rnd = gamrnd(a, b, method) % --*-- Unitary tests --*--

% This function produces independent random variates from the Gamma distribution.
%
%  INPUTS
%  - a       [double]    n*1 vector of positive parameters.
%  - b       [double]    n*1 vector of positive parameters.
%  - method  [struct]    Specifies which algorithms must be used.
%
%  OUTPUT
%  - rnd     [double]    n*1 vector of independent variates from the gamma(a,b) distribution.
%                        rnd(i) is gamma distributed with mean a(i)b(i) and variance a(i)b(i)^2.
%
%  REMARKS
%  The third input is a structure with two fields named `large` and `small`.
%  These fields define the algorithms to be used if a>1 (large) or a<1 (small).

% Copyright © 2006-2021 Dynare Team
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

%>
%>  Set defaults
%>  ------------

if nargin<2
    b = ones(size(a));
end

if nargin<3
    method = struct('large', 'Cheng', 'small', 'Johnk');
end

%>
%>  Check inputs
%>  ------------


[ma,na] = size(a);
[mb,nb] = size(b);

if ma~=mb || na~=nb
    error('gamrnd:: Input arguments must have the same size.');
end

if na~=1
    error('gamrnd:: Input arguments must be column vectors.');
end

if (any(a<0)) || (any(b<0)) || (any(a==Inf)) || (any(b==Inf))
    error('gamrnd:: Input arguments must be finite and positive.');
end

%>
%> Inititialize output
%> -------------------

rnd = NaN(ma,1);


% Get indices of integer (idx) and non integer (ddx) for the first hyperparameter a.
[~, idx, ddx] = isint(a);

if ~isempty(idx)
    % If the first hyperparameter (a) is an integer we can use the
    % exponential random number generator or rely in a Gaussian
    % approximation.
    sdx = find(a(idx)<30);
    ldx = find(a(idx)>=30);
    if ~isempty(sdx)
        % Exact sampling using random deviates from an exponential distribution.
        for i=1:length(sdx)
            rnd(idx(sdx(i))) = sum(exprnd(ones(a(idx(sdx(i))),1)))*b(idx(sdx(i)));
        end
    end
    if ~isempty(ldx)
        % Gaussian approximation.
        rnd(idx(ldx)) = sqrt(a(idx(ldx))).* b(idx(ldx)) .* randn(length(ldx), 1) + a(idx(ldx)) .* b(idx(ldx));
    end
end

if ~isempty(ddx)
    % The first hyperparameter is not an integer.
    sdx = find(a(ddx)<1);  % Indices for small a.
    ldx = find(a(ddx)>1);  % Indices for large a.
    if ~isempty(sdx)
        switch method.small
          case 'Weibull-rejection'
            rnd(ddx(sdx)) = gamrnd.weibull_rejection(a(ddx(sdx)),b(ddx(sdx)));
          case 'Johnk'
            rnd(ddx(sdx)) = gamrnd.johnk(a(ddx(sdx)),b(ddx(sdx)));
          case 'Berman'
            rnd(ddx(sdx)) = gamrnd.berman(a(ddx(sdx)),b(ddx(sdx)));
          case 'Ahrens-Dieter'
            rnd(ddx(sdx)) = gamrnd.ahrens_dieter(a(ddx(sdx)),b(ddx(sdx)));
          case 'Best'
            rnd(ddx(sdx)) = gamrnd.best_1983(a(ddx(sdx)),b(ddx(sdx)));
          otherwise
            error('Unknown algorithm for gamrnd.')
        end
    end
    if ~isempty(ldx)
        switch method.large
          case 'Knuth'
            rnd(ddx) = gamrnd.knuth(a(ddx),b(ddx));
          case 'Best'
            rnd(ddx(ldx)) = gamrnd.best_1978(a(ddx(ldx)),b(ddx(ldx)));
          case 'Cheng'
            rnd(ddx(ldx)) = gamrnd.cheng(a(ddx(ldx)),b(ddx(ldx)));
          otherwise
            error('Unknown algorithm for gamrnd.')
        end
    end
end

return

%@test:1
if ~isoctave && ~user_has_matlab_license('statistics_toolbox')
    method = struct('small', 'Weibull-rejection', 'large', 'Knuth');
    n = 1000000;
    m = 100;
    a = 0.1;
    b = 1.0;
    try
      mu = 0;
      s2 = 0;
      levels = .01:.01:10;
      ecdf = zeros(length(levels),1);
      for i = 1:m
          x = gamrnd(ones(n, 1)*a, ones(n,1)*b, method);
          mu = mu + mean(x);
          s2 = s2 + var(x);
          for j=1:length(levels)
              ecdf(j) = ecdf(j)+sum(x<levels(j))/n;
          end
      end
      mu = mu/m;
      s2 = s2/m;
      ecdf = ecdf/m;
      t(1) = true;
  catch
      t(1) = false;
  end
  if t(1)
      t(2) = abs(mu-a*b)<1e-3;
      t(3) = abs(s2-a*b^2)<1e-3;
      t(4) = max(abs(ecdf-gamcdf(transpose(levels), a, b)))<1e-3;
  end
  T = all(t);
else
    t = true(4, 1);
    T = true;
end
%@eof:1

%@test:2
if ~isoctave && ~user_has_matlab_license('statistics_toolbox')
    method = struct('small', 'Johnk', 'large', 'Knuth');
    n = 1000000;
    m = 100;
    a = 0.1;
    b = 1.0;
    try
        mu = 0;
        s2 = 0;
        levels = .01:.01:10;
        ecdf = zeros(length(levels),1);
        for i = 1:m
            x = gamrnd(ones(n, 1)*a, ones(n,1)*b, method);
            mu = mu + mean(x);
            s2 = s2 + var(x);
            for j=1:length(levels)
                ecdf(j) = ecdf(j)+sum(x<levels(j))/n;
            end
        end
        mu = mu/m;
        s2 = s2/m;
        ecdf = ecdf/m;
        t(1) = true;
    catch
        t(1) = false;
    end
    if t(1)
        t(2) = abs(mu-a*b)<1e-3;
        t(3) = abs(s2-a*b^2)<1e-3;
        t(4) = max(abs(ecdf-gamcdf(transpose(levels), a, b)))<1e-3;
    end
    T = all(t);
else
    t = true(4, 1);
    T = true;
end
%@eof:2

%@test:3
if ~isoctave && ~user_has_matlab_license('statistics_toolbox')
    method = struct('small', 'Berman', 'large', 'Knuth');
    n = 1000000;
    m = 100;
    a = 0.1;
    b = 1.0;
    try
        mu = 0;
        s2 = 0;
        levels = .01:.01:10;
        ecdf = zeros(length(levels),1);
        for i = 1:m
            x = gamrnd(ones(n, 1)*a, ones(n,1)*b, method);
            mu = mu + mean(x);
            s2 = s2 + var(x);
            for j=1:length(levels)
                ecdf(j) = ecdf(j)+sum(x<levels(j))/n;
            end
        end
        mu = mu/m;
        s2 = s2/m;
        ecdf = ecdf/m;
        t(1) = true;
    catch
        t(1) = false;
    end
    if t(1)
        t(2) = abs(mu-a*b)<1e-3;
        t(3) = abs(s2-a*b^2)<1e-3;
        t(4) = max(abs(ecdf-gamcdf(transpose(levels), a, b)))<1e-3;
    end
    T = all(t);
else
    t = true(4, 1);
    T = true;
end
%@eof:3

%@test:4
if ~isoctave && ~user_has_matlab_license('statistics_toolbox')
    method = struct('small', 'Ahrens-Dieter', 'large', 'Knuth');
    n = 1000000;
    m = 100;
    a = 0.1;
    b = 1.0;
    try
        mu = 0;
        s2 = 0;
        levels = .01:.01:10;
        ecdf = zeros(length(levels),1);
        for i = 1:m
            x = gamrnd(ones(n, 1)*a, ones(n,1)*b, method);
            mu = mu + mean(x);
            s2 = s2 + var(x);
            for j=1:length(levels)
                ecdf(j) = ecdf(j)+sum(x<levels(j))/n;
            end
        end
        mu = mu/m;
        s2 = s2/m;
        ecdf = ecdf/m;
        t(1) = true;
    catch
        t(1) = false;
    end
    if t(1)
        t(2) = abs(mu-a*b)<1e-3;
        t(3) = abs(s2-a*b^2)<1e-3;
        t(4) = max(abs(ecdf-gamcdf(transpose(levels), a, b)))<1e-3;
    end
    T = all(t);
else
    t = true(4, 1);
    T = true;
end
%@eof:4

%@test:5
if ~isoctave && ~user_has_matlab_license('statistics_toolbox')
    method = struct('small', 'Best', 'large', 'Knuth');
    n = 1000000;
    m = 100;
    a = 0.1;
    b = 1.0;
    try
        mu = 0;
        s2 = 0;
        levels = .01:.01:10;
        ecdf = zeros(length(levels),1);
        for i = 1:m
            x = gamrnd(ones(n, 1)*a, ones(n,1)*b, method);
            mu = mu + mean(x);
            s2 = s2 + var(x);
            for j=1:length(levels)
                ecdf(j) = ecdf(j)+sum(x<levels(j))/n;
            end
        end
        mu = mu/m;
        s2 = s2/m;
        ecdf = ecdf/m;
        t(1) = true;
    catch
        t(1) = false;
    end
    if t(1)
        t(2) = abs(mu-a*b)<1e-3;
        t(3) = abs(s2-a*b^2)<1e-3;
        t(4) = max(abs(ecdf-gamcdf(transpose(levels), a, b)))<1e-3;
    end
    T = all(t);
else
    t = true(4, 1);
    T = true;
end
%@eof:5

%@test:6
if ~isoctave && ~user_has_matlab_license('statistics_toolbox')
    method = struct('small', 'Weibull-rejection', 'large', 'Knuth');
    n = 1000000;
    m = 100;
    a = 1.5;
    b = 1.0;
    try
        mu = 0;
        s2 = 0;
        levels = .01:.01:15;
        ecdf = zeros(length(levels),1);
        for i = 1:m
            x = gamrnd(ones(n, 1)*a, ones(n,1)*b, method);
            mu = mu + mean(x);
            s2 = s2 + var(x);
            for j=1:length(levels)
                ecdf(j) = ecdf(j)+sum(x<levels(j))/n;
            end
        end
        mu = mu/m;
        s2 = s2/m;
        ecdf = ecdf/m;
        t(1) = true;
    catch
        t(1) = false;
    end
    if t(1)
        t(2) = abs(mu-a*b)<1e-3;
        t(3) = abs(s2-a*b^2)<1e-3;
        t(4) = max(abs(ecdf-gamcdf(transpose(levels), a, b)))<1e-3;
    end
    T = all(t);
else
    t = true(4, 1);
    T = true;
end
%@eof:6

%@test:7
if ~isoctave && ~user_has_matlab_license('statistics_toolbox')
    method = struct('small', 'Weibull-rejection', 'large', 'Cheng');
    n = 1000000;
    m = 100;
    a = 1.5;
    b = 1.0;
    try
        mu = 0;
        s2 = 0;
        levels = .01:.01:15;
        ecdf = zeros(length(levels),1);
        for i = 1:m
            x = gamrnd(ones(n, 1)*a, ones(n,1)*b, method);
            mu = mu + mean(x);
            s2 = s2 + var(x);
            for j=1:length(levels)
                ecdf(j) = ecdf(j)+sum(x<levels(j))/n;
            end
        end
        mu = mu/m;
        s2 = s2/m;
        ecdf = ecdf/m;
        t(1) = true;
    catch
        t(1) = false;
    end
    if t(1)
        t(2) = abs(mu-a*b)<1e-3;
        t(3) = abs(s2-a*b^2)<1e-3;
        t(4) = max(abs(ecdf-gamcdf(transpose(levels), a, b)))<1e-3;
    end
    T = all(t);
else
    t = true(4, 1);
    T = true;
end
%@eof:7

%@test:8
if ~isoctave && ~user_has_matlab_license('statistics_toolbox')
    method = struct('small', 'Weibull-rejection', 'large', 'Best');
    n = 1000000;
    m = 100;
    a = 1.5;
    b = 1.0;
    try
        mu = 0;
        s2 = 0;
        levels = .01:.01:15;
        ecdf = zeros(length(levels),1);
        for i = 1:m
            x = gamrnd(ones(n, 1)*a, ones(n,1)*b, method);
            mu = mu + mean(x);
            s2 = s2 + var(x);
            for j=1:length(levels)
                ecdf(j) = ecdf(j)+sum(x<levels(j))/n;
            end
        end
        mu = mu/m;
        s2 = s2/m;
        ecdf = ecdf/m;
        t(1) = true;
    catch
        t(1) = false;
    end
    if t(1)
        t(2) = abs(mu-a*b)<1e-3;
        t(3) = abs(s2-a*b^2)<1e-3;
        t(4) = max(abs(ecdf-gamcdf(transpose(levels), a, b)))<1e-3;
    end
    T = all(t);
else
    t = true(4, 1);
    T = true;
end
%@eof:8
