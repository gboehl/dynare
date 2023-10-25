function [y, T, success, max_res, iter] = solve_one_boundary(fh, y, x, params, steady_state, T, ...
                                                             y_index_eq, nze, periods, is_linear, Block_Num, y_kmin, maxit_, solve_tolf, cutoff, stack_solve_algo, is_forward, is_dynamic, verbose, M_, options_)
% Computes the deterministic simulation or the steady state for a block of equations containing
% only lags or only leads (but not both).
%
% INPUTS
%   fh                  [handle]        function handle to the static/dynamic file for the block
%   y                   [matrix]        All the endogenous variables of the model
%   x                   [matrix]        All the exogenous variables of the model
%   params              [vector]        All the parameters of the model
%   steady_state        [vector]        steady state of the model
%   T                   [matrix]        Temporary terms
%   y_index_eq          [vector of int] The index of the endogenous variables of
%                                       the block
%   nze                 [integer]       number of non-zero elements in the
%                                       jacobian matrix
%   periods             [integer]       number of simulation periods
%   is_linear           [logical]       whether the block is linear
%   Block_Num           [integer]       block number
%   y_kmin              [integer]       maximum number of lag in the model
%   maxit_              [integer]       maximum number of iteration in Newton
%   solve_tolf          [double]        convergence criteria
%   cutoff              [double]        cutoff to correct the direction in Newton in case
%                                       of singular jacobian matrix
%   stack_solve_algo    [integer]       linear solver method used in the Newton algorithm
%   is_forward          [logical]       Whether the block has to be solved forward
%                                       If false, the block is solved backward
%   is_dynamic          [logical]       Whether this is a deterministic simulation
%   verbose             [logical]       Whether iterations are to be printed
%   M_                  [structure]     storing the model information
%   options_            [structure]     storing the options
%
% OUTPUTS
%   y                   [matrix]        All endogenous variables of the model
%   T                   [matrix]        Temporary terms
%   success             [logical]       Whether a solution was found
%   max_res             [double]        ∞-norm of the residual
%   iter                [integer]       Number of iterations
%
% ALGORITHM
%   Newton with LU or GMRES or BicGstab for dynamic block

% Copyright © 1996-2023 Dynare Team
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

Blck_size=size(y_index_eq,2);
correcting_factor=0.01;
ilu_setup.type='crout';
ilu_setup.droptol=1e-10;
max_resa=1e100;
lambda = 1; % Length of Newton step
reduced = 0;
if is_forward
    incr = 1;
    start = y_kmin+1;
    finish = periods+y_kmin;
else
    incr = -1;
    start = periods+y_kmin;
    finish = y_kmin+1;
end

for it_=start:incr:finish
    cvg=false;
    iter=0;
    g1=spalloc( Blck_size, Blck_size, nze);
    while ~(cvg || iter>maxit_)
        if is_dynamic
            [yy, T(:, it_), r, g1] = fh(dynendo(y, it_, M_), x(it_, :), params, steady_state, ...
                                        M_.block_structure.block(Block_Num).g1_sparse_rowval, ...
                                        M_.block_structure.block(Block_Num).g1_sparse_colval, ...
                                        M_.block_structure.block(Block_Num).g1_sparse_colptr, T(:, it_));
            y(:, it_) = yy(M_.endo_nbr+(1:M_.endo_nbr));
        else
            [y, T, r, g1] = fh(y, x, params, M_.block_structure_stat.block(Block_Num).g1_sparse_rowval, ...
                               M_.block_structure_stat.block(Block_Num).g1_sparse_colval, ...
                               M_.block_structure_stat.block(Block_Num).g1_sparse_colptr, T);
        end
        if ~isreal(r)
            max_res=(-(max(max(abs(r))))^2)^0.5;
        else
            max_res=max(max(abs(r)));
        end
        if verbose
            disp(['iteration : ' int2str(iter+1) ' => ' num2str(max_res) ' time = ' int2str(it_)])
            if is_dynamic
                disp([char(M_.endo_names{y_index_eq}) repmat(' ', numel(y_index_eq), 1) num2str([y(y_index_eq, it_) r g1])])
            else
                disp([char(M_.endo_names{y_index_eq}) repmat(' ', numel(y_index_eq), 1) num2str([y(y_index_eq) r g1])])
            end
        end
        if ~isreal(max_res) || isnan(max_res)
            cvg = false;
        elseif is_linear && iter>0
            cvg = true;
        else
            cvg=(max_res<solve_tolf);
        end
        if ~cvg
            if iter>0
                if ~isreal(max_res) || isnan(max_res) || (max_resa<max_res && iter>1)
                    if isnan(max_res) || (max_resa<max_res && iter>0)
                        detJ=det(g1a);
                        if(abs(detJ)<1e-7)
                            max_factor=max(max(abs(g1a)));
                            ze_elem=sum(diag(g1a)<cutoff);
                            if verbose
                                disp([num2str(full(ze_elem),'%d') ' elements on the Jacobian diagonal are below the cutoff (' num2str(cutoff,'%f') ')'])
                            end
                            if correcting_factor<max_factor
                                correcting_factor=correcting_factor*4;
                                if verbose
                                    disp(['The Jacobian matrix is singular, det(Jacobian)=' num2str(detJ,'%f') '.'])
                                    disp(['    trying to correct the Jacobian matrix:'])
                                    disp(['    correcting_factor=' num2str(correcting_factor,'%f') ' max(Jacobian)=' num2str(full(max_factor),'%f')])
                                end
                                dx = - r/(g1+correcting_factor*speye(Blck_size));
                                y(y_index_eq, it_)=ya_save+lambda*dx;
                                continue
                            else
                                if verbose
                                    disp('The singularity of the jacobian matrix could not be corrected')
                                end
                                success = false;
                                return
                            end
                        end
                    elseif lambda>1e-8
                        lambda=lambda/2;
                        reduced = 1;
                        if verbose
                            disp(['reducing the path length: lambda=' num2str(lambda,'%f')])
                        end
                        if is_dynamic
                            y(y_index_eq, it_)=ya_save-lambda*dx;
                        else
                            y(y_index_eq)=ya_save-lambda*dx;
                        end
                        continue
                    else
                        if verbose
                            if cutoff==0
                                fprintf('Convergence not achieved in block %d, at time %d, after %d iterations.\n Increase "maxit".\n',Block_Num, it_, iter);
                            else
                                fprintf('Convergence not achieved in block %d, at time %d, after %d iterations.\n Increase "maxit" or set "cutoff=0" in model options.\n',Block_Num, it_, iter);
                            end
                        end
                        success = false;
                        return
                    end
                else
                    if lambda<1
                        lambda=max(lambda*2, 1);
                    end
                end
            end
            if is_dynamic
                ya = y(y_index_eq, it_);
            else
                ya = y(y_index_eq);
            end
            ya_save=ya;
            g1a=g1;
            if is_dynamic && stack_solve_algo==4
                stpmx = 100 ;
                stpmax = stpmx*max([sqrt(ya'*ya);size(y_index_eq,2)]);
                nn=1:size(y_index_eq,2);
                g = (r'*g1)';
                f = 0.5*r'*r;
                p = -g1\r ;
                [ya,f,r,check]=lnsrch1(ya,f,g,p,stpmax, ...
                                       @lnsrch1_wrapper_one_boundary,nn, ...
                                       nn, options_.solve_tolx, y_index_eq, fh, Block_Num, y, x, params, steady_state, T(:, it_), it_, M_);
                dx = ya - y(y_index_eq, it_);
                y(y_index_eq, it_) = ya;
                %% Recompute temporary terms, since they are not given as output of lnsrch1
                [~, T(:, it_)] = fh(dynendo(y, it_, M_), x(it_, :), params, steady_state, ...
                                    M_.block_structure.block(Block_Num).g1_sparse_rowval, ...
                                    M_.block_structure.block(Block_Num).g1_sparse_colval, ...
                                    M_.block_structure.block(Block_Num).g1_sparse_colptr, T(:, it_));
            elseif (is_dynamic && (stack_solve_algo==1 || stack_solve_algo==0 || stack_solve_algo==6)) || (~is_dynamic && options_.solve_algo==6)
                if verbose && ~is_dynamic
                    disp('steady: Sparse LU ')
                end
                dx =  g1\r;
                ya = ya - lambda*dx;
                if is_dynamic
                    y(y_index_eq, it_) = ya;
                else
                    y(y_index_eq) = ya;
                end
            elseif (stack_solve_algo==2 && is_dynamic) || (options_.solve_algo==7 && ~is_dynamic)
                flag1=1;
                if verbose && ~is_dynamic
                    disp('steady: GMRES ')
                end
                while flag1>0
                    [L1, U1]=ilu(g1,ilu_setup);
                    [dx,flag1] = gmres(g1,-r,Blck_size,1e-6,Blck_size,L1,U1);
                    if  flag1>0 || reduced
                        if verbose
                            if flag1==1
                                disp(['Error in simul: No convergence inside GMRES after ' num2str(iter,'%6d') ' iterations, in block' num2str(Block_Num,'%3d')])
                            elseif(flag1==2)
                                disp(['Error in simul: Preconditioner is ill-conditioned, in block' num2str(Block_Num,'%3d')])
                            elseif(flag1==3)
                                disp(['Error in simul: GMRES stagnated (Two consecutive iterates were the same.), in block' num2str(Block_Num,'%3d')])
                            end
                        end
                        ilu_setup.droptol = ilu_setup.droptol/10;
                        reduced = 0;
                    else
                        ya = ya + lambda*dx;
                        if is_dynamic
                            y(y_index_eq, it_) = ya;
                        else
                            y(y_index_eq) = ya';
                        end
                    end
                end
            elseif (stack_solve_algo==3 && is_dynamic) || (options_.solve_algo==8 && ~is_dynamic)
                flag1=1;
                if verbose && ~is_dynamic
                    disp('steady: BiCGStab')
                end
                while flag1>0
                    [L1, U1]=ilu(g1,ilu_setup);
                    [dx,flag1] = bicgstab(g1,-r,1e-6,Blck_size,L1,U1);
                    if flag1>0 || reduced
                        if verbose
                            if(flag1==1)
                                disp(['Error in simul: No convergence inside BiCGStab after ' num2str(iter,'%6d') ' iterations, in block' num2str(Block_Num,'%3d')])
                            elseif(flag1==2)
                                disp(['Error in simul: Preconditioner is ill-conditioned, in block' num2str(Block_Num,'%3d')])
                            elseif(flag1==3)
                                disp(['Error in simul: BiCGStab stagnated (Two consecutive iterates were the same.), in block' num2str(Block_Num,'%3d')])
                            end
                        end
                        ilu_setup.droptol = ilu_setup.droptol/10;
                        reduced = 0;
                    else
                        ya = ya + lambda*dx;
                        if is_dynamic
                            y(y_index_eq, it_) = ya;
                        else
                            y(y_index_eq) = ya';
                        end
                    end
                end
            else
                if is_dynamic
                    error(['options_.stack_solve_algo = ' num2str(stack_solve_algo) ' not implemented'])
                else
                    error(['options_.solve_algo = ' num2str(options_.solve_algo) ' not implemented'])
                end
            end
            iter=iter+1;
            max_resa = max_res;
        end
    end
    if ~cvg
        if verbose
            if cutoff == 0
                fprintf('Convergence not achieved in block %d, at time %d, after %d iterations.\n Increase "maxit".\n',Block_Num, it_,iter);
            else
                fprintf('Convergence not achieved in block %d, at time %d, after %d iterations.\n Increase "maxit" or set "cutoff=0" in model options.\n',Block_Num, it_,iter);
            end
        end
        success = false;
        return
    end
end

success = true;


function y3n = dynendo(y, it_, M_)
    if it_ > 1 && it_ < size(y, 2)
        y3n = reshape(y(:, it_+(-1:1)), 3*M_.endo_nbr, 1);
    elseif it_ > 1 % Purely backward model (in last period)
        y3n = [ reshape(y(:, it_+(-1:0)), 2*M_.endo_nbr, 1); NaN(M_.endo_nbr, 1) ];
    elseif it_ < size(y, 2) % Purely forward model (in first period)
        y3n = [ NaN(M_.endo_nbr, 1); reshape(y(:, it_+(0:1)), 2*M_.endo_nbr, 1) ];
    else % Static model
        y3n = [ NaN(M_.endo_nbr, 1); y(:, it_); NaN(M_.endo_nbr, 1) ]
    end

function r = lnsrch1_wrapper_one_boundary(ya, y_index, fh, Block_Num, y, x, params, steady_state, T, it_, M_)
    y(y_index, it_) = ya;
    [~, ~, r] = fh(dynendo(y, it_, M_), x(it_, :), params, steady_state, ...
                   M_.block_structure.block(Block_Num).g1_sparse_rowval, ...
                   M_.block_structure.block(Block_Num).g1_sparse_colval, ...
                   M_.block_structure.block(Block_Num).g1_sparse_colptr, T);
