/*
 * This file implements a small open economy with impatient agents, i.e. beta<1/(1+r),
 * but subject to a borrowing constraint d<2. The problem is set up as a 
 * mixed complementarity problem (MCP), requiring a setup of the form:
 *      LB =X     =>   F(X)>0,
 *      LB<=X<=UB =>   F(X)=0,
 *          X =UB =>   F(X)<0.
 * Thus, if d=2, the associated complementary Euler equation needs to return a 
 * negative residual. Economically, if the borrowing constraint binds, consumption 
 * today is lower than the unrestricted Euler equation 
 *
 *      1/c=beta*(1+r_bar)*1/c(+1)
 * 
 * would require and marginal utility on the left would be higher. Dynare by 
 * default computes the residual of an equation LHS=RHS as residual=LHS-RHS. For the 
 * above equation this would imply 
 *
 *      residual=1/c-beta*(1+r_bar)*1/c(+1)
 *   
 * and therefore would return a positive residual, while a negative one is required. 
 * The equation must therefore be entered as 
 *     beta*(1+r_bar)*1/c(+1)-1/c=0; 
 * which returns a negative residual when the upper bound is binding.
 *
 * This implementation was written by Johannes Pfeifer. If you spot any mistakes, 
 * email me at jpfeifer@gmx.de.
 * Please note that the following copyright notice only applies to this Dynare 
 * implementation of the model.
 */

/*
 * Copyright (C) 2022 Johannes Pfeifer
 *
 * This is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * It is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * For a copy of the GNU General Public License,
 * see <http://www.gnu.org/licenses/>.
 */

var  c d y;  
                                             
parameters  r_bar ${\bar r}$
            beta;
            
r_bar    = 0.04; %world interest rate		
beta=0.9;
d_bar=2;

model;
    [name='Evolution of debt']
    y + d= c + (1+r_bar)*d(-1);
    [name='Endowment']
    y = 1;
    [name='Euler equation']
    beta*(1+r_bar)*1/c(+1)-1/c=0 ⟂ d<2;
end;

initval;
d=d_bar;
y=1;
c=y-r_bar*d;
end;
resid;
steady(solve_algo=10);
resid;
