function y_=simult_(M_,options_,y0,dr,ex_,iorder)
% Simulates the model using a perturbation approach, given the path for the exogenous variables and the
% decision rules.
%
% INPUTS
%    M_       [struct]   model
%    options_ [struct]   options
%    y0       [double]   n*1 vector, initial value (n is the number of declared endogenous variables plus the number
%                        of auxilliary variables for lags and leads); must be in declaration order, i.e. as in M_.endo_names
%    dr       [struct]   matlab's structure where the reduced form solution of the model is stored.
%    ex_      [double]   T*q matrix of innovations.
%    iorder   [integer]  order of the taylor approximation.
%
% OUTPUTS
%    y_       [double]   n*(T+1) time series for the endogenous variables.
%
% SPECIAL REQUIREMENTS
%    none

% Copyright © 2001-2023 Dynare Team
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

iter = size(ex_,1);
endo_nbr = M_.endo_nbr;
exo_nbr = M_.exo_nbr;

if options_.loglinear && ~options_.logged_steady_state
    k = get_all_variables_but_lagged_leaded_exogenous(M_);
    dr.ys(k)=log(dr.ys(k));
end

if (iorder > 3) || (iorder == 3 && ~options_.pruning)
    if options_.order~=iorder
        error(['The k_order_simul routine requires the specified approximation order to be '...
               'consistent with the one used for computing the decision rules'])
    end
    y_ = k_order_simul(iorder,M_.nstatic,M_.npred,M_.nboth,M_.nfwrd,exo_nbr, ...
                       y0(dr.order_var,:),ex_',dr.ys(dr.order_var),dr, ...
                       options_.pruning);
    y_(dr.order_var,:) = y_;
else
    k2 = dr.kstate(find(dr.kstate(:,2) <= M_.maximum_lag+1),[1 2]);
    k2 = k2(:,1)+(M_.maximum_lag+1-k2(:,2))*endo_nbr;
    y_ = zeros(size(y0,1),iter+M_.maximum_lag);
    y_(:,1) = y0;
    order_var = dr.order_var;
    switch iorder
      case 1
        y_(:,1) = y_(:,1)-dr.ys;
        if isempty(dr.ghu)% For (linearized) deterministic models.
            for i = 2:iter+M_.maximum_lag
                yhat = y_(order_var(k2),i-1);
                y_(order_var,i) = dr.ghx*yhat;
            end
        elseif isempty(dr.ghx)% For (linearized) purely forward variables (no state variables).
            y_(dr.order_var,:) = dr.ghu*transpose(ex_);
        else
            epsilon = dr.ghu*transpose(ex_);
            for i = 2:iter+M_.maximum_lag
                yhat = y_(order_var(k2),i-1);
                y_(order_var,i) = dr.ghx*yhat + epsilon(:,i-1);
            end
        end
        y_ = bsxfun(@plus,y_,dr.ys);
      case 2
        constant = dr.ys(order_var)+.5*dr.ghs2;
        if options_.pruning
            y__ = y0;
            for i = 2:iter+M_.maximum_lag
                yhat1 = y__(order_var(k2))-dr.ys(order_var(k2));
                yhat2 = y_(order_var(k2),i-1)-dr.ys(order_var(k2));
                epsilon = ex_(i-1,:)';
                abcOut1 = A_times_B_kronecker_C(.5*dr.ghxx,yhat1);
                abcOut2 = A_times_B_kronecker_C(.5*dr.ghuu,epsilon);
                abcOut3 = A_times_B_kronecker_C(dr.ghxu,yhat1,epsilon);
                y_(order_var,i) = constant + dr.ghx*yhat2 + dr.ghu*epsilon ...
                    + abcOut1 + abcOut2 + abcOut3;
                y__(order_var) = dr.ys(order_var) + dr.ghx*yhat1 + dr.ghu*epsilon;
            end
        else
            for i = 2:iter+M_.maximum_lag
                yhat = y_(order_var(k2),i-1)-dr.ys(order_var(k2));
                epsilon = ex_(i-1,:)';
                abcOut1 = A_times_B_kronecker_C(.5*dr.ghxx,yhat);
                abcOut2 = A_times_B_kronecker_C(.5*dr.ghuu,epsilon);
                abcOut3 = A_times_B_kronecker_C(dr.ghxu,yhat,epsilon);
                y_(dr.order_var,i) = constant + dr.ghx*yhat + dr.ghu*epsilon ...
                    + abcOut1 + abcOut2 + abcOut3;
            end
        end
      case 3
        % only with pruning
        % the third moments of the shocks are assumed null. We don't have
        % an interface for specifying them
        ghx = dr.ghx;
        ghu = dr.ghu;
        ghxx = dr.ghxx;
        ghxu = dr.ghxu;
        ghuu = dr.ghuu;
        ghs2 = dr.ghs2;
        ghxxx = dr.ghxxx;
        ghxxu = dr.ghxxu;
        ghxuu = dr.ghxuu;
        ghuuu = dr.ghuuu;
        ghxss = dr.ghxss;
        ghuss = dr.ghuss;
        nspred = M_.nspred;
        ipred = M_.nstatic+(1:nspred);
        %construction follows Andreasen et al (2013), Technical
        %Appendix, Formulas (65) and (66)
        %split into first, second, and third order terms
        yhat1 = y0(order_var(k2))-dr.ys(order_var(k2));
        yhat2 = zeros(nspred,1);
        yhat3 = zeros(nspred,1);
        for i=2:iter+M_.maximum_lag
            u = ex_(i-1,:)';
            %construct terms of order 2 from second order part, based
            %on linear component yhat1
            gyy = A_times_B_kronecker_C(ghxx,yhat1);
            guu = A_times_B_kronecker_C(ghuu,u);
            gyu = A_times_B_kronecker_C(ghxu,yhat1,u);
            %construct terms of order 3 from second order part, based
            %on order 2 component yhat2
            gyy12 = A_times_B_kronecker_C(ghxx,yhat1,yhat2);
            gy2u = A_times_B_kronecker_C(ghxu,yhat2,u);
            %construct terms of order 3, all based on first order component yhat1
            y2a = kron(yhat1,yhat1);
            gyyy = A_times_B_kronecker_C(ghxxx,y2a,yhat1);
            u2a = kron(u,u);
            guuu = A_times_B_kronecker_C(ghuuu,u2a,u);
            yu = kron(yhat1,u);
            gyyu = A_times_B_kronecker_C(ghxxu,yhat1,yu);
            gyuu = A_times_B_kronecker_C(ghxuu,yu,u);
            %add all terms of order 3, linear component based on third
            %order yhat3
            yhat3 = ghx*yhat3 +gyy12 ... % prefactor is 1/2*2=1, see (65) Appendix Andreasen et al.
                    + gy2u ... % prefactor is 1/2*2=1, see (65) Appendix Andreasen et al.
                    + 1/6*(gyyy + guuu + 3*(gyyu + gyuu +  ghxss*yhat1 + ghuss*u)); %note: s is treated as variable, thus xss and uss are third order
            yhat2 = ghx*yhat2 + 1/2*(gyy + guu + 2*gyu + ghs2);
            yhat1 = ghx*yhat1 + ghu*u;
            y_(order_var,i) = dr.ys(order_var)+yhat1 + yhat2 + yhat3; %combine terms again
            yhat1 = yhat1(ipred);
            yhat2 = yhat2(ipred);
            yhat3 = yhat3(ipred);
        end
    end
end
