function N = RBC_MoM_steady_helper(THETA,ETAl,ETAc,BETTA,B,C_O_N,W)
if ETAc == 1 && ETAl == 1
    N = (1-BETTA*B)*(C_O_N*(1-B))^-1*W/THETA/(1+(1-BETTA*B)*(C_O_N*(1-B))^-1*W/THETA);
else
    % No closed-form solution use a fixed-point algorithm
    N0 = 1/3;
    N = fsolve(@(N) THETA*(1-N)^(-ETAl)*N^ETAc - (1-BETTA*B)*(C_O_N*(1-B))^(-ETAc)*W, N0,optimset('Display','off','TolX',1e-12,'TolFun',1e-12));
end