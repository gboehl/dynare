function [funcs, feedback_vars_idxs] = setup_time_recursive_block_simul(M_)
%function [funcs, feedback_vars_idxs] = setup_time_recursive_block_simul(M_)
%
% For solve_algo={12,14}, precompute the function handles for per-block dynamic files
% (it is surprisingly costly to recompute them within simulation loops).
% Also precompute indices of feedback variables (also brings some performance gains).
% By the way, do other sanity checks on block decomposition.

% Copyright © 2022 Dynare Team
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

if ~isfield(M_, 'block_structure')
    error('Can''t use solve_algo=12 nor solve_algo=14, because the block decomposition of the dynamic model failed')
end
if ~M_.block_structure.time_recursive
    error('Can''t use solve_algo=12 nor solve_algo=14, because the model is not purely backward/static/forward or you gave the ''block'' option to the ''model'' block')
end

funcs = cell(length(M_.block_structure.block), 1);
feedback_vars_idxs = cell(length(M_.block_structure.block), 1);
for blk = 1:length(M_.block_structure.block)
    funcs{blk} = str2func(sprintf('%s.sparse.block.dynamic_%d', M_.fname, blk));
    feedback_vars_idxs{blk} = M_.endo_nbr+M_.block_structure.block(blk).variable((M_.block_structure.block(blk).endo_nbr-M_.block_structure.block(blk).mfs+1):end); % Indices of feedback variables in the dynamic y vector (of size 3n)
end

end
