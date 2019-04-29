/*
 Check the policy functions obtained by perturbation at a high approximation
 order, using the Burnside (1998, JEDC) model (for which the analytical form of
 the policy function is known).

 As shown by Burnside, the policy function for yₜ is:

  yₜ =  βⁱ exp[aᵢ+bᵢ(xₜ−xₛₛ)]


 where:
                      θ²   ⎛     2ρ             1−ρ²ⁱ⎞
 — aᵢ = iθxₛₛ + σ² ─────── ⎢i − ────(1−ρⁱ) + ρ² ─────⎥
                   2(1−ρ)² ⎝    1−ρ             1−ρ² ⎠

         θρ
 — bᵢ = ───(1−ρⁱ)
        1−ρ

 — xₛₛ is the steady state of x
 — σ is the standard deviation of e.

 With some algebra, it can be shown that the derivative of yₜ at the deterministic
 steady state is equal to:

     ∂ᵐ⁺ⁿ⁺²ᵖ yₜ       ∞              (2p)!
  ──────────────── =  ∑  βⁱ bᵢᵐ⁺ⁿ ρᵐ ───── cᵢᵖ exp(iθxₛₛ)
  ∂ᵐxₜ₋₁ ∂ⁿeₜ ∂²ᵖs   ⁱ⁼¹               p!

 where:
 — s is the stochastic scale factor

           θ²   ⎛     2ρ             1−ρ²ⁱ⎞
 — cᵢ = ─────── ⎢i − ────(1−ρⁱ) + ρ² ─────⎥
        2(1−ρ)² ⎝    1−ρ             1−ρ² ⎠

 Note that derivatives with respect to an odd order for s (i.e. ∂²ᵖ⁺¹s) are always
 equal to zero.

 The policy function as returned in the oo_.dr.g_* matrices has the following properties:
 — its elements are pre-multiplied by the Taylor coefficients;
 — derivatives w.r.t. the stochastic scale factor have already been summed up;
 — symmetric elements are folded (and they are not pre-multiplied by the number of repetitions).

 As a consequence, the element gₘₙ corresponding to the m-th derivative w.r.t.
 to xₜ₋₁ and the n-th derivative w.r.t. to eₜ is given by:

          1                  ∞               cᵢᵖ
  gₘₙ = ──────      ∑        ∑  βⁱ bᵢᵐ⁺ⁿ ρᵐ ──── exp(iθxₛₛ)
        (m+n)!  0≤2p≤k-m-n  ⁱ⁼¹              p!

 where k is the order of approximation.

 */

@#define k = 9

var y x;
varexo e;

parameters beta theta rho xbar;
xbar = 0.0179;
rho = -0.139;
theta = -1.5;
theta = -10;
beta = 0.95;

model;
y = beta*exp(theta*x(+1))*(1+y(+1));
x = (1-rho)*xbar + rho*x(-1)+e;
end;

shocks;
var e; stderr 0.0348;
end;

initval;
x = xbar;
y = beta*exp(theta*xbar)/(1-beta*exp(theta*xbar));
end;

steady;

stoch_simul(order=@{k},k_order_solver,irf=0);

sigma2=M_.Sigma_e;
i = 1:800;
c = theta^2*sigma2/(2*(1-rho)^2)*(i-2*rho*(1-rho.^i)/(1-rho)+rho^2*(1-rho.^(2*i))/(1-rho^2));
b = theta*rho*(1-rho.^i)/(1-rho);

for ord = 0:@{k}
  g = oo_.dr.(['g_' num2str(ord)])(2,:); % Retrieve computed policy function for variable y

  for m = 0:ord % m is the derivation order with respect to x(-1)
    v = 0;
    for p = 0:floor((@{k}-ord)/2) % 2p is the derivation order with respect to s
      if ord+2*p > 0 % Skip the deterministic steady state constant
        v = v + sum(beta.^i.*exp(theta*xbar*i).*b.^ord.*rho^m.*c.^p)/factorial(ord)/factorial(p);
      end
    end
    if abs(v-g(ord+1-m)) > 1e-14
      error(['Error in matrix oo_.dr.g_' num2str(ord)])
    end
  end
end
