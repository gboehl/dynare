function [dr,info]=AIM_first_order_solver(g1,M_,dr,qz_criterium)

%@info:
%! @deftypefn {Function File} {[@var{dr},@var{info}] =} AIM_first_order_solver (@var{g1},@var{M_},@var{dr},@var{qz_criterium})
%! @anchor{AIM_first_order_solver}
%! @sp 1
%! Computes the first order reduced form of the DSGE model using AIM.
%! @sp 2
%! @strong{Inputs}
%! @sp 1
%! @table @ @var
%! @item g1
%! Matrix containing the Jacobian of the model
%! @item M_
%! Matlab's structure describing the model (initialized by @code{dynare}).
%! @item dr
%! Matlab's structure describing the reduced form solution of the model.
%! @item qz_criterium
%! Double containing the criterium to separate explosive from stable eigenvalues
%! @end table
%! @sp 2
%! @strong{Outputs}
%! @sp 1
%! @table @ @var
%! @item dr
%! Matlab's structure describing the reduced form solution of the model.
%! @item info
%! Integer scalar, error code.
%! @sp 1
%! @table @ @code
%! @item info==0
%! No error.
%! @item info==102
%! roots not correctly computed by real_schur
%! @item info==103
%! Blanchard & Kahn conditions are not satisfied: no stable equilibrium.
%! @item info==104
%! Blanchard & Kahn conditions are not satisfied: indeterminacy.
%! @item info==135
%! too many explosive roots and q(:,right) is singular
%! @item info==145
%! too few big roots, and q(:,right) is singular
%! @item info==105
%! q(:,right) is singular
%! @item info==161
%! too many exact siftrights
%! @item info==162
%! too many numeric shiftrights
%! @end table
%! @end table
%! @end deftypefn
%@eod:
%
% Maps Dynare jacobian to AIM 1st order model solver designed and developed by Gary Anderson
% and derives the solution for dr.ghx and dr.ghu from the AIM outputs
% AIM System is given as a sum:
% i.e. for i=-$...+&   SUM(Hi*xt+i)= £*zt, t = 0, . . . ,?
% and its input as single array of matrices: [H-$...  Hi ... H+&]
% and its solution as xt=SUM( Bi*xt+i) + @*£*zt for i=-$...-1
% with the output in form bb=[B-$...  Bi ... B-1] and @=inv(Ho+H1*B-1)
% Dynare jacobian = [fy'-$...  fy'i ... fy'+&  fu']
% where [fy'-$...  fy'i ... fy'+&]=[H-$...  Hi ... H+&] and fu'= £
%
% Dynare use:
%       1) the lognormal block in DR1 is being invoked for some models and changing
%       values of ghx and ghy. We need to return the AIM output
%       values before that block and run the block with the current returned values
%       of dr.ghx and dr.ghu if it is needed even when the AIM is used
%       (it does not depend on mjdgges output).
%
%       2) passing in aa={Q'|1}*g1 can produce ~ one order closer
%       results to the Dynare solutiion then when if plain g1 is passed,
%       i.e. diff < e-14 for aa and diff < *e-13 for g1 if Q' is used.
%
% Initially written by George Perendia

% Copyright © 2001-2024 Dynare Team
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

info = 0;

aimcode=-1;
neq = size(g1, 1); % no of equations
condn = 1.e-10; %SPAmalg uses this in zero tests
uprbnd = 1 + 1.e-6; %allow unit roots

try % try to run AIM
    [bb, rts, ~, ~, ~, ~, aimcode] = SPAmalg(g1(:, 1:3*M_.endo_nbr), neq, 1, 1, condn, uprbnd);
catch
    err = lasterror;
    disp(['Dynare AIM Solver error:' sprintf('%s; ID:%s',err.message, err.identifier)]);
    rethrow(lasterror);
end
if aimcode==1 %if OK
    dr.ghx = bb(dr.order_var, dr.order_var(M_.nstatic+(1:M_.nspred)));

    if M_.exo_nbr % if there are exogenous shocks then derive ghu for the shocks:
                  %   get H0 and H+1=HM
                  %    theH0= theAIM_H (:,M_.maximum_endo_lag*neq+1: (M_.maximum_endo_lag+1)*neq);
                  %theH0= theAIM_H (:,lags*neq+1: (lags+1)*neq);
                  %    theHP= theAIM_H (:,(M_.maximum_endo_lag+1)*neq+1: (M_.maximum_endo_lag+2)*neq);
                  %theHP= theAIM_H (:,(lags+1)*neq+1: (lags+2)*neq);
                                                                %? = inv(H0 + H1B1)
                                                                %phi= (theH0+theHP*sparse(bb(:,(lags-1)*neq+1:end)))\eye( neq);
                                                                %AIM_ghu=phi*theAIM_Psi;
                                                                %dr.ghu =AIM_ghu(dr.order_var,:); % order ghu
                                                                % Using AIM SPObstruct
        scof = SPObstruct(g1(:, 1:3*M_.endo_nbr), bb, neq, 1, 1);
        scof1 = scof(:, neq+1:end);
        scof1= scof1(:,dr.order_var);
        dr.ghu = -scof1 \ g1(:, 3*M_.endo_nbr+1:end);
    else
        dr.ghu = [];
    end
else
    err=SPAimerr(aimcode);
    disp(['Error in AIM: aimcode=' sprintf('%d : %s',aimcode, err)]);
    if aimcode < 1 || aimcode > 5  % too big exception, use mjdgges
        error('Error in AIM: aimcode=%d ; %s', aimcode, err);
    end
end

if aimcode ~=1
    info(1) = convertAimCodeToInfo(aimCode); %convert to be in the 100 range
    info(2) = 1.0e+8;
    return
end
A = kalman_transition_matrix(dr,M_.nstatic+(1:M_.nspred), 1:M_.nspred);
dr.eigval = eig(A);
disp(dr.eigval)
nd = M_.nspred+M_.nsfwrd;
nba = nd-sum( abs(dr.eigval) < qz_criterium );

nsfwrd = M_.nsfwrd;

if nba ~= nsfwrd
    temp = sort(abs(dr.eigval));
    if nba > nsfwrd
        temp = temp(nd-nba+1:nd-nsfwrd)-1-qz_criterium;
        info(1) = 3;
    elseif nba < nsfwrd
        temp = temp(nd-nsfwrd+1:nd-nba)-1-qz_criterium;
        info(1) = 4;
    end
    info(2) = temp'*temp;
    return
end


function [info] = convertAimCodeToInfo(aimCode)
% Returns an appropriate code for get_error_message.m (see that function for the meaning)

switch aimCode
  case 1
    info = 0;
  case 2
    info = 102;
  case 3
    info = 103;
  case 35
    info = 135;
  case 4
    info = 104;
  case 45
    info = 145;
  case 5
    info = 105;
  case 61
    info = 161;
  case 62
    info = 162;
  case 63
    info = 163;
  case 64
    info = 164;
  otherwise
    info = 1;
end
