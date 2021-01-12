function y = var_forecast(name, h, y, fcv)

% name : filename
% name        string    name of var model, provided in var statement
% h           int       number of steps-ahead forecast
% y           matrix    rows: realizations of endogenous variables in declaration order; cols: realizations in t, t-1, t-2 ... order of VAR
% fcv         string    name of variable we want forecast for

% returns the h-step-ahead VAR(order) forecast for fcv

% example calling:
% In Matlab:
% >> autoregressive_matrices{1} = [0.5000    0.1000; 0.4000    0.5000];
% >> autoregressive_matrices{2} = [0         0     ; 0.2500    0     ];
% >> mu                         = [0.0200; 0.0300];
% >> save('m1.mat', 'mu','autoregressive_matrices');

% In .mod file:
% var a b c d;
% ...
% var(model_name=m1,order=2) a c;

% From Matlab backend:
% >> yt   = [0.0600;    33.0000;    0.0300;    22.0000];
% >> ytm1 = [0.0550;    11.0000;    0.0300;    88.0000];
% >> var_forecast('m1', 1, [yt ytm1])
% >> var_forecast('m1', 2, [yt ytm1], ['a'])

%%
global M_;

%% construct y
assert( ...
    length(y) == length(M_.endo_names) ||             ... % when called from static model
    length(y) == sum(sum(M_.lead_lag_incidence ~= 0)) ... % when called from dynamic model
    );
yidx = zeros(size(M_.endo_names));
for i=1:size(M_.var.(name).var_list_,1)
    yidx = yidx | strcmp(strtrim(M_.var.(name).var_list_(i,:)), M_.endo_names);
end
y = y(yidx,:);

if nargin == 4
    fvidx = strcmp(fcv, M_.endo_names);
end

%% load .mat file
load(name, 'autoregressive_matrices', 'mu');
if ~exist('autoregressive_matrices', 'var') || ~exist('mu', 'var')
    error([name ' : must contain the variables autoregressive_matrices and mu']);
end
assert(h >= 1);

%% rewrite as VAR(1)
lm = length(mu);
lc = length(autoregressive_matrices);
assert(lc == M_.var.(name).order);
if size(y,1) ~= lm || size(y,2) ~= M_.var.(name).order
    error('The dimensions of y are not correct. It should be an nvars x order matrix');
end

A = zeros(lm*lc, lm*lc);
for i=1:lc
    if any([lm lm] ~= size(autoregressive_matrices{i}))
        error('The dimensions of mu and autoregressive_matrices are off');
    end
    col = lm*(i-1)+1:lm*i;
    A(1:lm, col) = autoregressive_matrices{i};
    if i ~= lc
        A(lm*i+1:lm*i+lm, col) = eye(lm, lm);
    end
end
if M_.var.(name).order > 1
    mu = [mu; zeros(lm*M_.var.(name).order-lm, 1)];
end

%% Calculate Forecast
%  New Introduction to Multiple Time Series Analysis
%  Helmut Lutkepohl
%  page 34
%
% An = eye(size(A));
% for i=1:h-1
%     An = An + A^i;
% end
% y = An*mu + A^h*y(:);

for i=1:h
    y = mu + A*y(:);
end
y = y(1:lm);

if nargin == 4
    retidx = find(fvidx & yidx == 1);
    if isempty(retidx)
        return;
    elseif retidx == 1
        y = y(1);
    else
        y = y(sum(yidx(1:retidx-1))+1);
    end
end
end