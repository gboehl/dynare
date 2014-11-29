function c = find_c(c0,pai,beta,epsilon,phi,gamma,omega)
c = csolve(@nk_ss,c0,[],1e-8,100,pai,beta,epsilon,phi,gamma,omega);
end

function r = nk_ss(c,pai,beta,epsilon,phi,gamma,omega)
r = pai*(pai-1)/c - beta*pai*(pai-1)/c-epsilon*phi*(c+(omega/2)*(pai-1)^2)^(gamma+1)/omega +(c+(omega/2)*(pai-1)^2)*(epsilon-1)/(omega*c);

end