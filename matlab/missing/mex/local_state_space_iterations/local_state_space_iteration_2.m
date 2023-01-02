function [y, y_] = local_state_space_iteration_2(yhat, epsilon, ghx, ghu, constant, ghxx, ghuu, ghxu, a, b, c)

% Given the demeaned states (yhat) and structural innovations (epsilon), this routine computes the level of selected endogenous variables when the
% model is approximated by an order two taylor expansion around the deterministic steady state. Depending on the number of input/output
% argument the pruning algorithm advocated by C. Sims is used or not (this version should not be used if the selected endogenous variables
% are not the states of the model).
%
% INPUTS
% - yhat        [double]    n×1 vector, initial condition for the state variables (centered).
% - epsilon     [double]    q×1 vector, innovations.
% - ghx         [double]    m×n matrix, restricted dr.ghx where we only consider the lines corresponding to a subset of endogenous variables.
% - ghu         [double]    m×q matrix, restricted dr.ghu where we only consider the lines corresponding to a subset of endogenous variables.
% - constant    [double]    m×1 vector, deterministic steady state plus second order correction for a subset of endogenous variables.
% - ghxx        [double]    m×n² matrix, restricted dr.ghxx where we only consider the lines corresponding to a subset of endogenous variables.
% - ghuu        [double]    m×q² matrix, restricted dr.ghuu where we only consider the lines corresponding to a subset of endogenous variables.
% - ghxu        [double]    m×nq matrix, restricted dr.ghxu where we only consider the lines corresponding to a subset of endogenous variables.
% - yhat_       [double]    [OPTIONAL] n×1 vector, initial condition for the state variables (centered) related to pruning.
% - ss          [double]    [OPTIONAL] n×1 vector, deterministic steady state for the state variables.
% - dummy       [double]    scalar, number of threads used.
%
% OUTPUTS
% - y           [double]    m×1 vector, (subset of the) endogenous variables.
% - y_          [double]    n×1 vector, update of the latent variables needed for the pruning version (first order update).
%
% REMARKS
% 1. If the function has more than nine input arguments (pruning) then it must have two output arguments (otherwise only one input).
% 2. Ninth input argument is not used, it is here only to have a common interface with the mex version.

% Copyright © 2011-2023 Dynare Team
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

if nargin==9
    pruning = false;
    if nargout>1
        error('local_state_space_iteration_2:: Numbers of input and output argument are inconsistent.')
    end
elseif nargin==11
    pruning = true;
    yhat_ = a;
    ss = b;
    if nargout~=2
        error('local_state_space_iteration_2:: Numbers of input and output argument are inconsistent.')
    end
else
    error('local_state_space_iteration_2:: Wrong number of input arguments.')
end

if pruning
    for i =1:size(yhat,2)
        y(:,i) = constant + ghx*yhat(:,i) + ghu*epsilon(:,i) ...
                 + A_times_B_kronecker_C(.5*ghxx, yhat_(:,i))  ...
                 + A_times_B_kronecker_C(.5*ghuu, epsilon(:,i)) ...
                 + A_times_B_kronecker_C(ghxu, yhat_(:,i), epsilon(:,i));
    end
    y_ = ghx*yhat_+ghu*epsilon;
    y_ = bsxfun(@plus, y_, ss);
else
    for i =1:size(yhat,2)
        y(:,i) = constant + ghx*yhat(:,i) + ghu*epsilon(:,i) ...
                 + A_times_B_kronecker_C(.5*ghxx, yhat(:,i))  ...
                 + A_times_B_kronecker_C(.5*ghuu, epsilon(:,i)) ...
                 + A_times_B_kronecker_C(ghxu, yhat(:,i), epsilon(:,i));
    end
end

return % --*-- Unit tests --*--

%@test:1
n = 2;
q = 3;

yhat = zeros(n,1);
epsilon = zeros(q,1);
ghx = rand(n,n);
ghu = rand(n,q);
constant = ones(n,1);
ghxx = rand(n,n*n);
ghuu = rand(n,q*q);
ghxu = rand(n,n*q);
yhat_ = zeros(n,1);
ss = ones(n,1);
addpath(sprintf('%s/missing/mex/local_state_space_iterations', fileparts(which('dynare'))))
% Call the tested routine.
for i=1:10
    y1 = local_state_space_iteration_2(yhat,epsilon,ghx,ghu,constant,ghxx,ghuu,ghxu,1);
    [y2,y2_] = local_state_space_iteration_2(yhat,epsilon,ghx,ghu,constant,ghxx,ghuu,ghxu,yhat_,ss,1);
end
rmpath(sprintf('%s/missing/mex/local_state_space_iterations', fileparts(which('dynare'))))
% Check the results.
t(1) = dassert(y1,ones(n,1));
t(2) = dassert(y2,ones(n,1));
t(3) = dassert(y2_,ones(n,1));
T = all(t);
%@eof:1

%@test:2
old_path = pwd;
cd([fileparts(which('dynare')) '/../tests/particle']);
load dsgebase2data;
cd(old_path);
n = length(state_variables_idx);
m = rows(ReducedForm.ghx);
q = size(ReducedForm.ghu,2);
yhatinit = randn(n,1);
epsilon = randn(q,2);
ghx = ReducedForm.ghx;
ghu = ReducedForm.ghu;
constant = ReducedForm.constant;
ghxx = ReducedForm.ghxx;
ghuu = ReducedForm.ghuu;
ghxu = ReducedForm.ghxu;
yhatinit_ = randn(n,1);
steadystate = ReducedForm.steadystate;
t = true(6,1);
% Call the tested routine (matlab).
addpath(sprintf('%s/missing/mex/local_state_space_iterations', fileparts(which('dynare'))))
try
    yhat1 = local_state_space_iteration_2(yhatinit, epsilon(:,1), ghx, ghu, constant, ghxx, ghuu, ghxu, 1);
    yhat1 = local_state_space_iteration_2(yhat1(state_variables_idx), epsilon(:,2), ghx, ghu, constant, ghxx, ghuu, ghxu, 1);
catch
    t(1) = false;
end
try
    [yhat2, yhat2_] = local_state_space_iteration_2(yhatinit, epsilon(:,1), ghx, ghu, constant, ghxx, ghuu, ghxu, yhatinit_, steadystate, 1);
    [yhat2, yhat2_] = local_state_space_iteration_2(yhat2(state_variables_idx), epsilon(:,2), ghx, ghu, constant, ghxx, ghuu, ghxu, yhat2_(state_variables_idx), steadystate, 1);
catch
    t(2) = false;
end
rmpath(sprintf('%s/missing/mex/local_state_space_iterations', fileparts(which('dynare'))))
% Call the tested routine (mex).
try
    yhat3 = local_state_space_iteration_2(yhatinit, epsilon(:,1), ghx, ghu, constant, ghxx, ghuu, ghxu, 1);
    yhat3 = local_state_space_iteration_2(yhat3(state_variables_idx), epsilon(:,2), ghx, ghu, constant, ghxx, ghuu, ghxu, 1);
catch
    t(3) = false;
end
try
    [yhat4, yhat4_] = local_state_space_iteration_2(yhatinit, epsilon(:,1), ghx, ghu, constant, ghxx, ghuu, ghxu, yhatinit_, steadystate, 1);
    [yhat4, yhat4_] = local_state_space_iteration_2(yhat4(state_variables_idx), epsilon(:,2), ghx, ghu, constant, ghxx, ghuu, ghxu, yhat4_(state_variables_idx), steadystate, 1);
catch
    t(4) = false;
end
t(5) = max(abs(yhat1-yhat3))<1e-12; % Compare matlab and mex routines without pruning.
t(6) = max(abs(yhat2-yhat4))<1e-12; % Compare matlab and mex routines with pruning.
                                    % Check the results.
T = all(t);
%@eof:2

%@test:3
if false
    % TIMING TEST (parallelization with openmp)
    old_path = pwd;
    cd([fileparts(which('dynare')) '/../tests/particle']);
    load dsgebase2data;
    cd(old_path);
    n = length(state_variables_idx);
    q = size(ReducedForm.ghu,2);
    yhat = randn(n,10000000);
    epsilon = .01*randn(q,10000000);
    ghx = ReducedForm.ghx(state_variables_idx,:);
    ghu = ReducedForm.ghu(state_variables_idx,:);
    constant = ReducedForm.constant(state_variables_idx,:);
    ghxx = ReducedForm.ghxx(state_variables_idx,:);
    ghuu = ReducedForm.ghuu(state_variables_idx,:);
    ghxu = ReducedForm.ghxu(state_variables_idx,:);
    yhatinit_ = randn(n,1);
    ss = ReducedForm.state_variables_steady_state;
    yhat_ = randn(n,10000000);
    t = NaN(4,1);
    tic, for i=1:10, y1 = local_state_space_iteration_2(yhat,epsilon,ghx,ghu,constant,ghxx,ghuu,ghxu,1); end
    t1 = toc;
    tic, for i=1:10, y2 = local_state_space_iteration_2(yhat,epsilon,ghx,ghu,constant,ghxx,ghuu,ghxu,2); end
    t2 = toc;
    t(1) = dassert(y1,y2,1e-15); clear('y1');
    tic, for i=1:10, y3 = local_state_space_iteration_2(yhat,epsilon,ghx,ghu,constant,ghxx,ghuu,ghxu,3); end
    t3 = toc;
    t(2) = dassert(y2,y3,1e-15); clear('y2');
    tic, for i=1:10, y4 = local_state_space_iteration_2(yhat,epsilon,ghx,ghu,constant,ghxx,ghuu,ghxu,4); end
    t4 = toc;
    t(3) = dassert(y4,y3,1e-15); clear('y3','y4');
    t(4) = (t1>t2) && (t2>t3) && (t3>t4);
    if ~t(4)
        disp('Timmings:')
        [t1, t2, t3, t4]
    end
    % Check the results.
    T = all(t);
else
    t = 1;
    T = true;
end
%@eof:3