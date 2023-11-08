function [y, T, success, err, iter] = solve_two_boundaries_lbj(fh, y, x, steady_state, T, blk, options_, M_)
% Computes the deterministic simulation of a block of equations containing
% both lead and lag variables, using the LBJ algorithm.
%
% INPUTS
%   fh                  [handle]        function handle to the dynamic file for the block
%   y                   [matrix]        All the endogenous variables of the model
%   x                   [matrix]        All the exogenous variables of the model
%   steady_state        [vector]        steady state of the model
%   T                   [matrix]        Temporary terms
%   blk                 [integer]       block number
%   options_            [structure]     storing the options
%   M_                  [structure]     Model description
%
% OUTPUTS
%   y                   [matrix]        All endogenous variables of the model
%   T                   [matrix]        Temporary terms
%   success             [logical]       Whether a solution was found
%   err                 [double]        ∞-norm of Δy
%   iter                [integer]       Number of iterations
%
% ALGORITHM
%   Laffargue, Boucekkine, Juillard (LBJ)
%   see Juillard (1996) Dynare: A program for the resolution and
%   simulation of dynamic models with forward variables through the use
%   of a relaxation algorithm. CEPREMAP. Couverture Orange. 9602.

% Copyright © 2023 Dynare Team
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

sparse_rowval = M_.block_structure.block(blk).g1_sparse_rowval;
sparse_colval = M_.block_structure.block(blk).g1_sparse_colval;
sparse_colptr = M_.block_structure.block(blk).g1_sparse_colptr;

periods = options_.periods;

% NB: notations are deliberately similar to those of sim1_lbj.m

ny = M_.block_structure.block(blk).mfs;

% Compute which columns, in the 3×n-wide Jacobian, have non-zero elements
% corresponding to the forward (iyf) or backward (iyp) variables
iyp = find(sparse_colptr(2:ny+1)-sparse_colptr(1:ny));
iyf = find(sparse_colptr(2*ny+2:end)-sparse_colptr(2*ny+1:end-1));

y_index = M_.block_structure.block(blk).variable(end-ny+1:end);

success = false;

for iter = 1:options_.simul.maxit
    h = clock;

    c = zeros(ny*periods, length(iyf)+1); % Stores the D and d of Sébastien’s presentation
    it_ = M_.maximum_lag+1;
    [yy, T(:, it_), r, g1] = fh(dynendo(y, it_, M_), x(it_, :), M_.params, steady_state, ...
                                sparse_rowval, sparse_colval, sparse_colptr, T(:, it_));
    y(:, it_) = yy(M_.endo_nbr+(1:M_.endo_nbr));
    ic = 1:ny;
    icp = iyp;
    c(ic, :) = full(g1(:, ny+(1:ny))) \ [ full(g1(:, 2*ny+iyf)) -r ];
    for it_ = M_.maximum_lag+(2:periods)
        [yy, T(:, it_), r, g1] = fh(dynendo(y, it_, M_), x(it_, :), M_.params, steady_state, ...
                                    sparse_rowval, sparse_colval, sparse_colptr, T(:, it_));
        y(:, it_) = yy(M_.endo_nbr+(1:M_.endo_nbr));
        j = [ full(g1(:, ny+(1:ny))) -r ];
        j(:, [ iyf ny+1 ]) = j(:, [ iyf ny+1 ]) - full(g1(:, iyp)) * c(icp, :);
        ic = ic + ny;
        icp = icp + ny;
        c(ic, :) = j(:, 1:ny) \ [ full(g1(:, 2*ny+iyf)) j(:, ny+1) ];
    end

    dy = back_subst_lbj(c, ny, iyf, periods);

    y(y_index, M_.maximum_lag+(1:periods)) = y(y_index, M_.maximum_lag+(1:periods)) + dy;
    err = max(max(abs(dy)));

    if options_.verbosity
        fprintf('Iter: %s,\t err. = %s, \t time = %s\n', num2str(iter), num2str(err), num2str(etime(clock, h)));
    end
    if err < options_.dynatol.x
        success = true;
        break
    end
end


function y3n = dynendo(y, it_, M_)
    y3n = reshape(y(:, it_+(-1:1)), 3*M_.endo_nbr, 1);
