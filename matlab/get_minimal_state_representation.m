function [CheckCO,minnx,minA,minB,minC,minD,dminA,dminB,dminC,dminD] = get_minimal_state_representation(A,B,C,D,dA,dB,dC,dD)
% [CheckCO,minnx,minA,minB,minC,minD,dminA,dminB,dminC,dminD] = get_minimal_state_representation(A,B,C,D,dA,dB,dC,dD)
% Derives and checks the minimal state representation of the ABCD
% representation of a state space model
% -------------------------------------------------------------------------
% INPUTS
%   A:          [endo_nbr by endo_nbr] Transition matrix from Kalman filter
%                 for all endogenous declared variables, in DR order
%   B:          [endo_nbr by exo_nbr]  Transition matrix from Kalman filter
%                 mapping shocks today to endogenous variables today, in DR order
%   C:          [obs_nbr by endo_nbr] Measurement matrix from Kalman filter
%                 linking control/observable variables to states, in DR order
%   D:          [obs_nbr by exo_nbr] Measurement matrix from Kalman filter
%                 mapping shocks today to controls/observables today, in DR order
%   dA:         [endo_nbr by endo_nbr by totparam_nbr] in DR order
%                 Jacobian (wrt to all parameters) of transition matrix A
%   dB:         [endo_nbr by exo_nbr by totparam_nbr] in DR order
%                 Jacobian (wrt to all parameters) of transition matrix B
%   dC:         [obs_nbr by endo_nbr by totparam_nbr] in DR order
%                 Jacobian (wrt to all parameters) of measurement matrix C
%   dD:         [obs_nbr by exo_nbr by totparam_nbr] in DR order
%                 Jacobian (wrt to all parameters) of measurement matrix D
% -------------------------------------------------------------------------
% OUTPUTS
%   CheckCO:    [scalar] indicator, equals to 1 if minimal state representation is found
%   minnx:      [scalar] length of minimal state vector
%   minA:       [minnx by minnx] Transition matrix A for evolution of minimal state vector
%   minB:       [minnx by exo_nbr] Transition matrix B for evolution of minimal state vector
%   minC:       [obs_nbr by minnx] Measurement matrix C for evolution of controls, depending on minimal state vector only
%   minD:       [obs_nbr by minnx] Measurement matrix D for evolution of controls, depending on minimal state vector only
%   dminA:      [minnx by minnx by totparam_nbr] in DR order
%                 Jacobian (wrt to all parameters) of transition matrix minA
%   dminB:      [minnx by exo_nbr by totparam_nbr] in DR order
%                 Jacobian (wrt to all parameters) of transition matrix minB
%   dminC:      [obs_nbr by minnx by totparam_nbr] in DR order
%                 Jacobian (wrt to all parameters) of measurement matrix minC
%   dminD:      [obs_nbr by exo_nbr by totparam_nbr] in DR order
%                 Jacobian (wrt to all parameters) of measurement matrix minD
% -------------------------------------------------------------------------
% This function is called by 
%   * get_identification_jacobians.m (previously getJJ.m)
% -------------------------------------------------------------------------
% This function calls
%   * check_minimality (embedded)
% =========================================================================
% Copyright (C) 2019 Dynare Team
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
% =========================================================================

nx = size(A,2);
ny = size(C,1);
nu = size(B,2);

% Check controllability and observability conditions for full state vector
CheckCO = check_minimality(A,B,C);

if CheckCO == 1 % If model is already minimal    
    minnx = nx;
    minA = A; 
    minB = B; 
    minC = C; 
    minD = D;
    if nargout > 6
        dminA = dA;
        dminB = dB;
        dminC = dC;
        dminD = dD;
    end
else
    %Model is not minimal
    realsmall = 1e-7;
    % create indices for unnecessary states
    endogstateindex = find(abs(sum(A,1))<realsmall);
    exogstateindex = find(abs(sum(A,1))>realsmall);
    minnx = length(exogstateindex);    
    % remove unnecessary states from solution matrices
    A_1 = A(endogstateindex,exogstateindex);
    A_2 = A(exogstateindex,exogstateindex);
    B_2 = B(exogstateindex,:); 
    C_1 = C(:,endogstateindex); 
    C_2 = C(:,exogstateindex);
    D   = D;
    ind_A1 = endogstateindex;
    ind_A2 = exogstateindex;
    % minimal representation
    minA = A_2; 
    minB = B_2;
    minC = C_2;
    minD = D;        
    % Check controllability and observability conditions
    CheckCO = check_minimality(minA,minB,minC);

    if CheckCO ~=1
        j=1;
        while (CheckCO==0 && j<minnx)
            combis = nchoosek(1:minnx,j);
            i=1;
            while i <= size(combis,1)
                ind_A2 = exogstateindex;
                ind_A1 = [endogstateindex ind_A2(combis(j,:))];
                ind_A2(combis(j,:)) = [];
                % remove unnecessary states from solution matrices
                A_1 = A(ind_A1,ind_A2);
                A_2 = A(ind_A2,ind_A2);
                B_2 = B(ind_A2,:); 
                C_1 = C(:,ind_A1); 
                C_2 = C(:,ind_A2);
                D = D;
                % minimal representation
                minA = A_2; 
                minB = B_2;
                minC = C_2;
                minD = D;
                % Check controllability and observability conditions
                CheckCO = check_minimality(minA,minB,minC);
                if CheckCO == 1                    
                    minnx = length(ind_A2);
                    break;                    
                end
                i=i+1;
            end
            j=j+1;
        end
    end
    if nargout > 6
        dminA = dA(ind_A2,ind_A2,:);
        dminB = dB(ind_A2,:,:);
        dminC = dC(:,ind_A2,:);
        dminD = dD;
    end
end

function CheckCO = check_minimality(A,B,C)
%% This function computes the controllability and the observability matrices of the ABCD system and checks if the system is minimal
%
% Inputs: Solution matrices A,B,C of ABCD representation of a DSGE model
% Outputs: CheckCO: equals 1, if both rank conditions for observability and controllability are fullfilled
n = size(A,1);
CC = B; % Initialize controllability matrix
OO = C; % Initialize observability matrix
if n >= 2
    for indn = 1:1:n-1
        CC = [CC, (A^indn)*B]; % Set up controllability matrix
        OO = [OO; C*(A^indn)]; % Set up observability matrix
    end
end
CheckC = (rank(full(CC))==n);   % Check rank of controllability matrix
CheckO = (rank(full(OO))==n);   % Check rank of observability matrix
CheckCO = CheckC&CheckO;        % equals 1 if minimal state
end % check_minimality end

end % Main function end
