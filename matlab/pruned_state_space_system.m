function dr = pruned_state_space_system(M, options, dr)
% This can be set up into the following ABCD representation:
% z = c + A*z(-1) + B*xi
% y = ybar + d + C*z(-1) + D*xi

% Denote
%   hx  = dr.ghx(indx,:);  hu  = dr.ghu(indx,:);
%   hxx = dr.ghxx(indx,:); hxu = dr.ghxu(indx,:); huu = dr.ghuu(indx,:); hs2 = dr.ghs2(indx,:);
%   gx  = dr.ghx(indy,:);    gu  = dr.ghu(indy,:);
%   gxx = dr.ghxx(indy,:);   gxu = dr.ghxu(indy,:);   guu = dr.ghuu(indy,:);   gs2 = dr.ghs2(indy,:);
%   hxss = dr.ghxss(indx,:); huss = dr.ghuss(indx,:); hxxx = dr.ghxxx(indx,:); huuu = dr.ghuuu(indx,:); hxxu = dr.ghxxu(indx,:); hxuu = dr.ghxxu(indx,:);
%   gxss = dr.ghxss(indy,:);   guss = dr.ghuss(indy,:);   gxxx = dr.ghxxx(indy,:);   guuu = dr.ghuuu(indy,:);   gxxu = dr.ghxxu(indy,:);   gxuu = dr.ghxxu(indy,:);
% Law of motion for first-order effects states xf and selected endogenous first-order effects variables yf: 
%   xf = hx*xf(-1) + hu*u
%   yf = gx*xf(-1) + gu*u
% Law of motion for second-order effects states xs and selected endogenous second-order effects variables ys:
%   xs = hx*xs(-1) + 1/2*hxx*kron(xf(-1),xf(-1)) + 1/2*huu*kron(u,u) + hxu*kron(xf(-1),u) + 1/2*hs2
%   ys = gx*xs(-1) + 1/2*gxx*kron(xf(-1),xf(-1)) + 1/2*guu*kron(u,u) + gxu*kron(xf(-1),u) + 1/2*gs2
% Law of motion for third-order effects states xrd and selected endogenous second-order effects variables yrd:
%   xrd = hx*xrd(-1) + hxx*kron(xf(-1),xs(-1)) + hxu*kron(xs(-1),u) + 3/6*hxss*xf(-1) + 3/6*huss*u + 1/6*hxxx*kron(xf(-1),kron(xf(-1),xf(-1))) + 1/6*huuu*kron(u,kron(u,u)) + 3/6*hxxu*kron(xf(-1),kron(xf(-1),u)) + 3/6*hxuu*kron(xf(-1),kron(u,u))
%   yrd = gx*xrd(-1) + gxx*kron(xf(-1),xs(-1)) + gxu*kron(xs(-1),u) + 3/6*gxss*xf(-1) + 3/6*guss*u + 1/6*gxxx*kron(xf(-1),kron(xf(-1),xf(-1))) + 1/6*guuu*kron(u,kron(u,u)) + 3/6*gxxu*kron(xf(-1),kron(xf(-1),u)) + 3/6*gxuu*kron(xf(-1),kron(u,u))
    % Selected variables:       y  = ybar + yf + ys
    % Pruned state vector:      z  = [xf; xs; kron(xf,xf)]
    % Pruned innovation vector: xi = [u; kron(u,u)-vec(Sigma_u); kron(xf(-1),u); kron(u,xf(-1))]    
    % Second moments of xi
    % Sig_xi = E[u*u'        ,               u*kron(u,u)'                           ,   zeros(u_nbr,x_nbr*u_nbr),          zeros(u_nbr,u_nbr*x_nbr);
    %            kron(u,u)*u',               kron(u,u)*kron(u,u)'-Sig_u(:)*Sig_u(:)',   zeros(u_nbr^2,x_nbr*u_nbr),        zeros(u_nbr^2,u_nbr*x_nbr);
    %            zeros(x_nbr*u_nbr,u_nbr),   zeros(x_nbr*u_nbr,u_nbr^2),                kron(xf(-1),u)*kron(xf(-1),u)',    kron(xf(-1),u)*kron(u,xf(-1))';
    %            zeros(u_nbr*x_nbr,u_nbr),   zeros(u_nbr*x_nbr,u_nbr^2),                kron(u,xf(-1))*kron(xf(-1),u)',    kron(u,xf(-1))*kron(u,xf(-1))';]
    % That is, we only need to compute:
    % - Sig_xi_11 = E[u*u']=Sig_u
    % - Sig_xi_21 = E[kron(u,u)*u']=0 for symmetric distributions
    % - Sig_xi_12 = transpose(Sig_xi_21)= 0 for symmetric distributions
    % - Sig_xi_22 = E[kron(u,u)*kron(u,u)']-Sig_u(:)*Sig_u(:)' which requires the fourth-order product moments of u
    % - Sig_xi_34 = E[kron(xf(-1),u)*kron(xf(-1),u)'] which requires the second-order moments of xf and second-order moments of u
    % - Sig_xi_43 = transpose(Sig_xi_34)
    % - Sig_xi_33 = 
    % - Sig_xi_44 = 


%% Auxiliary indices, numbers and flags
order = options.order;
indx = [M.nstatic+(1:M.nspred) M.endo_nbr+(1:size(dr.ghx,2)-M.nspred)]';
indy = (1:M.endo_nbr)'; %by default select all variables
compute_derivs = false;
if isfield(options,'options_ident')
    compute_derivs  = true;
    stderrparam_nbr = length(options.options_ident.indpstderr);
    corrparam_nbr   = size(options.options_ident.indpcorr,1);
    modparam_nbr    = length(options.options_ident.indpmodel);
    totparam_nbr    = stderrparam_nbr+corrparam_nbr+modparam_nbr;
    indy            = options.options_ident.indvobs;
end
Yss = dr.ys(dr.order_var);                  if compute_derivs; dYss = dr.derivs.dYss; end
E_uu = M.Sigma_e;                           if compute_derivs; dE_uu = dr.derivs.dSigma_e; end
Correlation_matrix = M.Correlation_matrix;  if compute_derivs; dCorrelation_matrix = dr.derivs.dCorrelation_matrix; end
ghx = dr.ghx;                               if compute_derivs; dghx = dr.derivs.dghx; end
ghu = dr.ghu;                               if compute_derivs; dghu = dr.derivs.dghu; end
Om  = ghu*E_uu*transpose(ghu);           if compute_derivs; dOm = dr.derivs.dOm; end
y_nbr = length(indy);
x_nbr = length(indx);
u_nbr = M.exo_nbr;

% indices for extended state vector z and extended shock vector e
id_z1_xf = (1:x_nbr);
id_e1_u  = (1:u_nbr);
if options.order > 1
    id_z2_xs    = id_z1_xf(end)   + (1:x_nbr);
    id_z3_xf_xf = id_z2_xs(end)   + (1:x_nbr*x_nbr);    
    id_e2_u_u   = id_e1_u(end)    + (1:u_nbr*u_nbr);
    id_e3_xf_u  = id_e2_u_u(end)  + (1:x_nbr*u_nbr);
    id_e4_u_xf  = id_e3_xf_u(end) + (1:u_nbr*x_nbr);
    
    ghxx = dr.ghxx;  if compute_derivs; dghxx = dr.derivs.dghxx; end
    ghxu = dr.ghxu;  if compute_derivs; dghxu = dr.derivs.dghxu; end
    ghuu = dr.ghuu;  if compute_derivs; dghuu = dr.derivs.dghuu; end
    ghs2 = dr.ghs2;  if compute_derivs; dghs2 = dr.derivs.dghs2; end
end
if options.order > 2
    id_z4_xrd      = id_z3_xf_xf(end)   + (1:x_nbr);
    id_z5_xf_xs    = id_z4_xrd(end)     + (1:x_nbr*x_nbr);
    id_z6_xf_xf_xf = id_z5_xf_xs(end)   + (1:x_nbr*x_nbr*x_nbr);    
    id_e5_xs_u     = id_e4_u_xf(end)    + (1:x_nbr*u_nbr);
    id_e6_u_xs     = id_e5_xs_u(end)    + (1:u_nbr*x_nbr);
    id_e7_xf_xf_u  = id_e6_u_xs(end)    + (1:x_nbr*x_nbr*u_nbr);
    id_e8_xf_u_xf  = id_e7_xf_xf_u(end) + (1:x_nbr*u_nbr*x_nbr);
    id_e9_u_xf_xf  = id_e8_xf_u_xf(end) + (1:u_nbr*x_nbr*x_nbr);
    id_e10_xf_u_u  = id_e9_u_xf_xf(end) + (1:x_nbr*u_nbr*u_nbr);
    id_e11_u_xf_u  = id_e10_xf_u_u(end) + (1:u_nbr*x_nbr*u_nbr);
    id_e12_u_u_xf  = id_e11_u_xf_u(end) + (1:u_nbr*u_nbr*x_nbr);
    id_e13_u_u_u   = id_e12_u_u_xf(end) + (1:u_nbr*u_nbr*u_nbr);
    
    ghxxx = dr.ghxxx; if compute_derivs; dghxxx = dr.derivs.dghxxx; end
    ghxxu = dr.ghxxu; if compute_derivs; dghxxu = dr.derivs.dghxxu; end
    ghxuu = dr.ghxuu; if compute_derivs; dghxuu = dr.derivs.dghxuu; end
    ghuuu = dr.ghuuu; if compute_derivs; dghuuu = dr.derivs.dghuuu; end
    ghxss = dr.ghxss; if compute_derivs; dghxss = dr.derivs.dghxss; end
    ghuss = dr.ghuss; if compute_derivs; dghuss = dr.derivs.dghuss; end
end

%% First-order
z_nbr = x_nbr;
e_nbr = M.exo_nbr;
A = ghx(indx,:);
B = ghu(indx,:);
C = ghx(indy,:);
D = ghu(indy,:);
c = zeros(x_nbr,1);
d = zeros(y_nbr,1);
Varinov = E_uu;
E_z = zeros(x_nbr,1); %this is E_xf
if compute_derivs == 1
    dA = dghx(indx,:,:);
    dB = dghu(indx,:,:);
    dC = dghx(indy,:,:);
    dD = dghu(indy,:,:);
    dc = zeros(x_nbr,totparam_nbr);
    dd = zeros(y_nbr,totparam_nbr);
    dVarinov = dE_uu;
    dE_z = zeros(x_nbr,totparam_nbr);
end

if order > 1
    % Some common and useful objects for order > 1
    % Compute E[xf*xf']
    E_xfxf = lyapunov_symm(ghx(indx,:), Om(indx,indx), options.lyapunov_fixed_point_tol, options.qz_criterium, options.lyapunov_complex_threshold, 1, options.debug); %use 1 to initialize persistent variables
    K_ux = commutation(u_nbr,x_nbr);    %commutation matrix
    K_xu = commutation(x_nbr,u_nbr);    %commutation matrix
    QPu  = quadruplication(u_nbr);      %quadruplication matrix as in Meijer (2005), similar to Magnus-Neudecker definition of duplication matrix (i.e. dyn_vec) but for fourth-order product moments
    %Compute unique product moments of E[kron(u*u',uu')] = E[kron(u,u)*kron(u,u)']
    COMBOS4 = flipud(allVL1(u_nbr, 4)); %all possible combinations of powers that sum up to four for fourth-order product moments of u
    E_u_u_u_u = zeros(u_nbr*(u_nbr+1)/2*(u_nbr+2)/3*(u_nbr+3)/4,1); %only unique entries of E[kron(u,kron(u,kron(u,u)))]
    if compute_derivs && (stderrparam_nbr+corrparam_nbr>0)
        dE_u_u_u_u = zeros(u_nbr*(u_nbr+1)/2*(u_nbr+2)/3*(u_nbr+3)/4,stderrparam_nbr+corrparam_nbr);
    end
    for j = 1:size(COMBOS4,1)
        E_u_u_u_u(j) = prodmom(E_uu, 1:u_nbr, COMBOS4(j,:));
        if compute_derivs && (stderrparam_nbr+corrparam_nbr>0)
            dE_u_u_u_u(j,1:(stderrparam_nbr+corrparam_nbr)) = prodmom_deriv(E_uu, 1:u_nbr, COMBOS4(j,:), dE_uu(:,:,1:(stderrparam_nbr+corrparam_nbr)), dCorrelation_matrix(:,:,1:(stderrparam_nbr+corrparam_nbr)));
        end
    end
    E_u_u_u_u = QPu*E_u_u_u_u; %add duplicate product moments, i.e. this is E[kron(u,kron(u,kron(u,u)))]
    E_uu_uu = reshape(E_u_u_u_u,u_nbr^2,u_nbr^2); %E[kron(u*u',uu')] = E[kron(u,u)*kron(u,u)']

    % E[kron((xf*xf'),(u*u'))]
    E_xfxf_E_uu = kron(E_xfxf,E_uu);
    if compute_derivs
        dE_xfxf = zeros(size(E_xfxf,1),size(E_xfxf,2),totparam_nbr);
        dE_xfxf_E_uu = zeros(size(E_xfxf_E_uu,1),size(E_xfxf_E_uu,2),totparam_nbr);
        for jp = 1:totparam_nbr
            if jp <= (stderrparam_nbr+corrparam_nbr)
                %Jacobian of E_xfxf wrt exogenous paramters: dE_xfxf(:,:,jp)-hx*dE_xfxf(:,:,jp)*hx'=dOm(:,:,jp), because dhx is zero by construction for stderr and corr parameters
                dE_xfxf(:,:,jp) = lyapunov_symm(ghx(indx,:), dOm(indx,indx,jp),options.lyapunov_fixed_point_tol,options.qz_criterium,options.lyapunov_complex_threshold,2,options.debug);
            else 
                %Jacobian of E_xfxf wrt model parameters: dE_xfxf(:,:,jp) - hx*dE_xfxf(:,:,jp)*hx' = d(hu*Sig_u*hu')(:,:,jp) + dhx(:,:,jp)*E_xfxf*hx'+ hx*E_xfxf*dhx(:,:,jp)'
                dE_xfxf(:,:,jp) = lyapunov_symm(ghx(indx,:), dghx(indx,:,jp)*E_xfxf*ghx(indx,:)'+ghx(indx,:)*E_xfxf*dghx(indx,:,jp)'+dOm(indx,indx,jp),options.lyapunov_fixed_point_tol,options.qz_criterium,options.lyapunov_complex_threshold,2,options.debug);
                %here method=2 is used to spare a lot of computing time while not repeating Schur every time                
            end            
            dE_xfxf_E_uu(:,:,jp) = kron(dE_xfxf(:,:,jp),E_uu) + kron(E_xfxf,dE_uu(:,:,jp));
        end        
    end
    
    % Second-order pruned state space system
    z_nbr = x_nbr+x_nbr+x_nbr*x_nbr;
    e_nbr = u_nbr+u_nbr*u_nbr+x_nbr*u_nbr+u_nbr*x_nbr;
    
    A = zeros(z_nbr, z_nbr);
    A(id_z1_xf    , id_z1_xf   ) = ghx(indx,:);
    A(id_z2_xs    , id_z2_xs   ) = ghx(indx,:);
    A(id_z2_xs    , id_z3_xf_xf) = 0.5*ghxx(indx,:);
    A(id_z3_xf_xf , id_z3_xf_xf) = kron(ghx(indx,:),ghx(indx,:));
    
    B = zeros(z_nbr, e_nbr);
    B(id_z1_xf    , id_e1_u   ) = ghu(indx,:);
    B(id_z2_xs    , id_e2_u_u ) = 0.5*ghuu(indx,:);
    B(id_z2_xs    , id_e3_xf_u) = ghxu(indx,:);
    B(id_z3_xf_xf , id_e2_u_u ) = kron(ghu(indx,:),ghu(indx,:));
    B(id_z3_xf_xf , id_e3_xf_u) = kron(ghx(indx,:),ghu(indx,:));
    B(id_z3_xf_xf , id_e4_u_xf) = kron(ghu(indx,:),ghx(indx,:));    
    
    C = zeros(y_nbr, z_nbr);
    C(1:y_nbr , id_z1_xf   ) = ghx(indy,:);
    C(1:y_nbr , id_z2_xs   ) = ghx(indy,:);
    C(1:y_nbr , id_z3_xf_xf) = 0.5*ghxx(indy,:);    
    
    D = zeros(y_nbr, e_nbr);
    D(1:y_nbr , id_e1_u   ) = ghu(indy,:);
    D(1:y_nbr , id_e2_u_u ) = 0.5*ghuu(indy,:);
    D(1:y_nbr , id_e3_xf_u) = ghxu(indy,:);
    
    c = zeros(z_nbr, 1);
    c(id_z2_xs    , 1) = 0.5*ghs2(indx,:) + 0.5*ghuu(indx,:)*E_uu(:);
    c(id_z3_xf_xf , 1) = kron(ghu(indx,:),ghu(indx,:))*E_uu(:);
    
    d = zeros(y_nbr, 1);
    d(1:y_nbr, 1) = 0.5*ghs2(indy,:) + 0.5*ghuu(indy,:)*E_uu(:);

    Varinov = zeros(e_nbr,e_nbr);
    Varinov(id_e1_u    , id_e1_u)    = E_uu;
    Varinov(id_e2_u_u  , id_e2_u_u)  = E_uu_uu-E_uu(:)*transpose(E_uu(:));
    Varinov(id_e3_xf_u , id_e3_xf_u) = E_xfxf_E_uu;
    Varinov(id_e4_u_xf , id_e3_xf_u) = K_ux*E_xfxf_E_uu;
    Varinov(id_e3_xf_u , id_e4_u_xf) = transpose(Varinov(id_e4_u_xf,id_e3_xf_u));
    Varinov(id_e4_u_xf , id_e4_u_xf) = K_ux*E_xfxf_E_uu*transpose(K_ux);
    
    I_hx    = speye(x_nbr)-ghx(indx,:);
    I_hxinv = I_hx\speye(x_nbr);
    E_xs = I_hxinv*(0.5*ghxx(indx,:)*E_xfxf(:) + c(x_nbr+(1:x_nbr),1));
    
    E_z = [zeros(x_nbr,1); E_xs; E_xfxf(:)];
    if compute_derivs        
        dA = zeros(size(A,1),size(A,2),totparam_nbr);
        dB = zeros(size(B,1),size(B,2),totparam_nbr);
        dC = zeros(size(C,1),size(C,2),totparam_nbr);
        dD = zeros(size(D,1),size(D,2),totparam_nbr);
        dc = zeros(size(c,1),totparam_nbr);
        dd = zeros(size(d,1),totparam_nbr);
        dVarinov = zeros(size(Varinov,1),size(Varinov,2),totparam_nbr);
        dE_xs = zeros(size(E_xs,1),size(E_xs,2),totparam_nbr);
        dE_z = zeros(size(E_z,1),size(E_z,2),totparam_nbr);
        for jp = 1:totparam_nbr
            dA(id_z1_xf    , id_z1_xf    , jp) = dghx(indx,:,jp);
            dA(id_z2_xs    , id_z2_xs    , jp) = dghx(indx,:,jp);
            dA(id_z2_xs    , id_z3_xf_xf , jp) = 0.5*dghxx(indx,:,jp);
            dA(id_z3_xf_xf , id_z3_xf_xf , jp) = kron(dghx(indx,:,jp),ghx(indx,:)) + kron(ghx(indx,:),dghx(indx,:,jp));

            dB(id_z1_xf    , id_e1_u    , jp) = dghu(indx,:,jp);
            dB(id_z2_xs    , id_e2_u_u  , jp) = 0.5*dghuu(indx,:,jp);
            dB(id_z2_xs    , id_e3_xf_u , jp) = dghxu(indx,:,jp);
            dB(id_z3_xf_xf , id_e2_u_u  , jp) = kron(dghu(indx,:,jp),ghu(indx,:)) + kron(ghu(indx,:),dghu(indx,:,jp));
            dB(id_z3_xf_xf , id_e3_xf_u , jp) = kron(dghx(indx,:,jp),ghu(indx,:)) + kron(ghx(indx,:),dghu(indx,:,jp));
            dB(id_z3_xf_xf , id_e4_u_xf , jp) = kron(dghu(indx,:,jp),ghx(indx,:)) + kron(ghu(indx,:),dghx(indx,:,jp));

            dC(1:y_nbr , id_z1_xf    , jp) = dghx(indy,:,jp);
            dC(1:y_nbr , id_z2_xs    , jp) = dghx(indy,:,jp);
            dC(1:y_nbr , id_z3_xf_xf , jp) = 0.5*dghxx(indy,:,jp);

            dD(1:y_nbr , id_e1_u    , jp) = dghu(indy,:,jp);
            dD(1:y_nbr , id_e2_u_u  , jp) = 0.5*dghuu(indy,:,jp);
            dD(1:y_nbr , id_e3_xf_u , jp) = dghxu(indy,:,jp);

            dc(id_z2_xs    , jp) = 0.5*dghs2(indx,jp) + 0.5*(dghuu(indx,:,jp)*E_uu(:) + ghuu(indx,:)*vec(dE_uu(:,:,jp)));
            dc(id_z3_xf_xf , jp) = (kron(dghu(indx,:,jp),ghu(indx,:)) + kron(ghu(indx,:),dghu(indx,:,jp)))*E_uu(:) + kron(ghu(indx,:),ghu(indx,:))*vec(dE_uu(:,:,jp));

            dd(1:y_nbr , jp) = 0.5*dghs2(indy,jp) + 0.5*(dghuu(indy,:,jp)*E_uu(:) + ghuu(indy,:)*vec(dE_uu(:,:,jp)));

            if jp <= (stderrparam_nbr+corrparam_nbr)
                dVarinov(id_e1_u   , id_e1_u   , jp) = dE_uu(:,:,jp);
                dVarinov(id_e2_u_u , id_e2_u_u , jp) = reshape(QPu*dE_u_u_u_u(:,jp), u_nbr^2, u_nbr^2) - (reshape(dE_uu(:,:,jp),u_nbr^2,1)*transpose(E_uu(:)) + E_uu(:)*transpose(reshape(dE_uu(:,:,jp),u_nbr^2,1)));
            end
            dVarinov(id_e3_xf_u , id_e3_xf_u , jp) = dE_xfxf_E_uu(:,:,jp);
            dVarinov(id_e4_u_xf , id_e3_xf_u , jp) = K_ux*dE_xfxf_E_uu(:,:,jp);
            dVarinov(id_e3_xf_u , id_e4_u_xf , jp) = transpose(dVarinov(id_e4_u_xf,id_e3_xf_u,jp));
            dVarinov(id_e4_u_xf , id_e4_u_xf , jp) = dVarinov(id_e4_u_xf , id_e3_xf_u , jp)*transpose(K_ux);
            
            dE_xs(:,jp) = I_hxinv*( dghx(indx,:,jp)*E_xs + 1/2*(dghxx(indx,:,jp)*E_xfxf(:) + ghxx(indx,:)*vec(dE_xfxf(:,:,jp)) + dc(x_nbr+(1:x_nbr),jp)) );
            dE_z(:,jp) = [zeros(x_nbr,1); dE_xs(:,jp); vec(dE_xfxf(:,:,jp))];
        end
    end

    if order > 2
        Q6Pu  = Q6_plication(u_nbr);      %quadruplication matrix as in Meijer (2005), similar to Magnus-Neudecker definition of duplication matrix (i.e. dyn_vec) but for fourth-order product moments
        %Compute unique product moments of E[kron(u*u',uu')] = E[kron(u,u)*kron(u,u)']
        COMBOS6 = flipud(allVL1(u_nbr, 6)); %all possible combinations of powers that sum up to four for fourth-order product moments of u
        E_u_u_u_u_u_u = zeros(u_nbr*(u_nbr+1)/2*(u_nbr+2)/3*(u_nbr+3)/4*(u_nbr+4)/5*(u_nbr+5)/6,1); %only unique entries of E[kron(u,kron(u,kron(u,kron(u,kron(u,u)))))]
        if compute_derivs && (stderrparam_nbr+corrparam_nbr>0)
            dE_u_u_u_u_u_u = zeros(size(E_u_u_u_u_u_u,1),stderrparam_nbr+corrparam_nbr);
        end
        for j = 1:size(COMBOS6,1)
            E_u_u_u_u_u_u(j) = prodmom(E_uu, 1:u_nbr, COMBOS6(j,:));
            if compute_derivs && (stderrparam_nbr+corrparam_nbr>0)
                dE_u_u_u_u_u_u(j,1:(stderrparam_nbr+corrparam_nbr)) = prodmom_deriv(E_uu, 1:u_nbr, COMBOS6(j,:), dE_uu(:,:,1:(stderrparam_nbr+corrparam_nbr)), dCorrelation_matrix(:,:,1:(stderrparam_nbr+corrparam_nbr)));
            end
        end

        Var_z = lyapunov_symm(A, B*Varinov*transpose(B), options.lyapunov_fixed_point_tol, options.qz_criterium, options.lyapunov_complex_threshold, 1, options.debug);
        %E_xfxf = Var_z(id_z1_xf,id_z1_xf); %this is E_xfxf, as E_xf=0, we already have that
        E_xfxs = Var_z(id_z1_xf,id_z2_xs); %as E_xf=0
        E_xf_xf_xf = vec(Var_z(id_z1_xf,id_z3_xf_xf)); %as E_xf=0, this is vec(E[xf*kron(xf,xf)'])
        E_xsxf = Var_z(id_z2_xs,id_z1_xf); %as E_xf=0        
        E_xsxs = Var_z(id_z2_xs,id_z2_xs) + E_xs*transpose(E_xs);
        E_xs_xf_xf = vec(Var_z(id_z2_xs,id_z3_xf_xf)+E_xs*E_xfxf(:)'); %this is vec(E[xs*kron(xf,xf)'])
        %E_xf_xf_xf = vec(Var_z(id_z3_xf_xf,id_z1_xf));
        E_xf_xf_xs = vec(Var_z(id_z3_xf_xf,id_z2_xs) + E_xfxf(:)*E_xs');
        E_xf_xf_xf_xf = vec(Var_z(id_z3_xf_xf,id_z3_xf_xf) + E_xfxf(:)*E_xfxf(:)');
        %E_xs_E_uu = kron(E_xs,E_uu);        
        E_xfxs_E_uu = kron(E_xfxs,E_uu);
        
        z_nbr = x_nbr+x_nbr+x_nbr^2+x_nbr+x_nbr*x_nbr+x_nbr^3;
        e_nbr = u_nbr+u_nbr^2+x_nbr*u_nbr+u_nbr*x_nbr+x_nbr*u_nbr+u_nbr*x_nbr+x_nbr*x_nbr*u_nbr+x_nbr*u_nbr*x_nbr+u_nbr*x_nbr*x_nbr+u_nbr*u_nbr*x_nbr+u_nbr*x_nbr*u_nbr+x_nbr*u_nbr*u_nbr+u_nbr*u_nbr*u_nbr;
        hx_hx = kron(ghx(indx,:),ghx(indx,:));
        hu_hu = kron(ghu(indx,:),ghu(indx,:));
        hx_hu = kron(ghx(indx,:),ghu(indx,:));
        hu_hx = kron(ghu(indx,:),ghx(indx,:));
        
        if compute_derivs
            dVar_z = zeros(size(Var_z,1),size(Var_z,2),totparam_nbr);
            dE_xfxs = zeros(size(E_xfxs,1),size(E_xfxs,2),totparam_nbr);
            dE_xf_xf_xf = zeros(size(E_xf_xf_xf,1),totparam_nbr);
            dE_xsxf = zeros(size(E_xsxf,1),size(E_xsxf,2),totparam_nbr);
            dE_xsxs = zeros(size(E_xsxs,1),size(E_xsxs,2),totparam_nbr);
            dE_xs_xf_xf = zeros(size(E_xs_xf_xf,1),totparam_nbr);
            dE_xf_xf_xs = zeros(size(E_xf_xf_xs,1),totparam_nbr);
            dE_xf_xf_xf_xf = zeros(size(E_xf_xf_xf_xf,1),totparam_nbr);
            dE_xfxs_E_uu = zeros(size(E_xfxs_E_uu,1),size(E_xfxs_E_uu,2),totparam_nbr);
            dhx_hx = zeros(size(hx_hx,1),size(hx_hx,2),totparam_nbr);
            dhu_hu = zeros(size(hu_hu,1),size(hu_hu,2),totparam_nbr);
            dhx_hu = zeros(size(hx_hu,1),size(hx_hu,2),totparam_nbr);
            dhu_hx = zeros(size(hu_hx,1),size(hu_hx,2),totparam_nbr);
            for jp=1:totparam_nbr
                dVar_z(:,:,jp) = lyapunov_symm(A, dB(:,:,jp)*Varinov*transpose(B) + B*dVarinov(:,:,jp)*transpose(B) +B*Varinov*transpose(dB(:,:,jp)) + dA(:,:,jp)*Var_z*A' + A*Var_z*transpose(dA(:,:,jp)),... ,...
                                               options.lyapunov_fixed_point_tol, options.qz_criterium, options.lyapunov_complex_threshold, 2, options.debug); %2 is used as we use 1 above
                dE_xfxs(:,:,jp) = dVar_z(id_z1_xf,id_z2_xs,jp);
                dE_xf_xf_xf(:,jp) = vec(dVar_z(id_z1_xf,id_z3_xf_xf,jp));
                dE_xsxf(:,:,jp) = dVar_z(id_z2_xs,id_z1_xf,jp);
                dE_xsxs(:,:,jp) = dVar_z(id_z2_xs,id_z2_xs,jp) + dE_xs(:,jp)*transpose(E_xs) + E_xs*transpose(dE_xs(:,jp));
                dE_xs_xf_xf(:,jp) = vec(dVar_z(id_z2_xs,id_z3_xf_xf,jp) + dE_xs(:,jp)*E_xfxf(:)' + E_xs*vec(dE_xfxf(:,:,jp))');
                dE_xf_xf_xs(:,jp) = vec(dVar_z(id_z3_xf_xf,id_z2_xs,jp) + vec(dE_xfxf(:,:,jp))*E_xs' + E_xfxf(:)*dE_xs(:,jp)');
                dE_xf_xf_xf_xf(:,jp) = vec(dVar_z(id_z3_xf_xf,id_z3_xf_xf,jp) + vec(dE_xfxf(:,:,jp))*E_xfxf(:)' + E_xfxf(:)*vec(dE_xfxf(:,:,jp))');
                dE_xfxs_E_uu(:,:,jp) = kron(dE_xfxs(:,:,jp),E_uu) + kron(E_xfxs,dE_uu(:,:,jp));

                dhx_hx(:,:,jp) = kron(dghx(indx,:,jp),ghx(indx,:)) + kron(ghx(indx,:),dghx(indx,:,jp));
                dhu_hu(:,:,jp) = kron(dghu(indx,:,jp),ghu(indx,:)) + kron(ghu(indx,:),dghu(indx,:,jp));
                dhx_hu(:,:,jp) = kron(dghx(indx,:,jp),ghu(indx,:)) + kron(ghx(indx,:),dghu(indx,:,jp));
                dhu_hx(:,:,jp) = kron(dghu(indx,:,jp),ghx(indx,:)) + kron(ghu(indx,:),dghx(indx,:,jp));
            end
        end

        A = zeros(z_nbr,z_nbr);
        A(id_z1_xf       , id_z1_xf      ) = ghx(indx,:);
        A(id_z2_xs       , id_z2_xs      ) = ghx(indx,:);
        A(id_z2_xs       , id_z3_xf_xf   ) = 1/2*ghxx(indx,:);
        A(id_z3_xf_xf    , id_z3_xf_xf   ) = hx_hx;
        A(id_z4_xrd      , id_z1_xf      ) = 3/6*ghxss(indx,:);
        A(id_z4_xrd      , id_z4_xrd     ) = ghx(indx,:);
        A(id_z4_xrd      , id_z5_xf_xs   ) = ghxx(indx,:);
        A(id_z4_xrd      , id_z6_xf_xf_xf) = 1/6*ghxxx(indx,:);
        A(id_z5_xf_xs    , id_z1_xf      ) = kron(ghx(indx,:),1/2*ghs2(indx,:));
        A(id_z5_xf_xs    , id_z5_xf_xs   ) = hx_hx;
        A(id_z5_xf_xs    , id_z6_xf_xf_xf) = kron(ghx(indx,:),1/2*ghxx(indx,:));
        A(id_z6_xf_xf_xf , id_z6_xf_xf_xf) = kron(ghx(indx,:),hx_hx);

        B = zeros(z_nbr,e_nbr);
        B(id_z1_xf       , id_e1_u      ) = ghu(indx,:);
        B(id_z2_xs       , id_e2_u_u    ) = 1/2*ghuu(indx,:);
        B(id_z2_xs       , id_e3_xf_u   ) = ghxu(indx,:);
        B(id_z3_xf_xf    , id_e2_u_u    ) = hu_hu;
        B(id_z3_xf_xf    , id_e3_xf_u   ) = hx_hu;
        B(id_z3_xf_xf    , id_e4_u_xf   ) = hu_hx;
        B(id_z4_xrd      , id_e1_u      ) = 3/6*ghuss(indx,:);
        B(id_z4_xrd      , id_e5_xs_u   ) = ghxu(indx,:);
        B(id_z4_xrd      , id_e7_xf_xf_u) = 3/6*ghxxu(indx,:);
        B(id_z4_xrd      , id_e10_xf_u_u) = 3/6*ghxuu(indx,:);
        B(id_z4_xrd      , id_e13_u_u_u ) = 1/6*ghuuu(indx,:);
        B(id_z5_xf_xs    , id_e1_u      ) = kron(ghu(indx,:),1/2*ghs2(indx,:));
        B(id_z5_xf_xs    , id_e6_u_xs   ) = hu_hx;
        B(id_z5_xf_xs    , id_e7_xf_xf_u) = kron(ghx(indx,:),ghxu(indx,:));
        B(id_z5_xf_xs    , id_e9_u_xf_xf) = kron(ghu(indx,:),1/2*ghxx(indx,:));
        B(id_z5_xf_xs    , id_e10_xf_u_u) = kron(ghx(indx,:),1/2*ghuu(indx,:));
        B(id_z5_xf_xs    , id_e11_u_xf_u) = kron(ghu(indx,:),ghxu(indx,:));
        B(id_z5_xf_xs    , id_e13_u_u_u ) = kron(ghu(indx,:),1/2*ghuu(indx,:));
        B(id_z6_xf_xf_xf , id_e7_xf_xf_u) = kron(hx_hx,ghu(indx,:));
        B(id_z6_xf_xf_xf , id_e8_xf_u_xf) = kron(ghx(indx,:),hu_hx);
        B(id_z6_xf_xf_xf , id_e9_u_xf_xf) = kron(ghu(indx,:),hx_hx);
        B(id_z6_xf_xf_xf , id_e10_xf_u_u) = kron(hx_hu,ghu(indx,:));
        B(id_z6_xf_xf_xf , id_e11_u_xf_u) = kron(ghu(indx,:),hx_hu);
        B(id_z6_xf_xf_xf , id_e12_u_u_xf) = kron(hu_hu,ghx(indx,:));
        B(id_z6_xf_xf_xf , id_e13_u_u_u ) = kron(ghu(indx,:),hu_hu);

        C = zeros(y_nbr,z_nbr);
        C(1:y_nbr , id_z1_xf      ) = ghx(indy,:) + 1/2*ghxss(indy,:);
        C(1:y_nbr , id_z2_xs      ) = ghx(indy,:);
        C(1:y_nbr , id_z3_xf_xf   ) = 1/2*ghxx(indy,:);
        C(1:y_nbr , id_z4_xrd     ) = ghx(indy,:);
        C(1:y_nbr , id_z5_xf_xs   ) = ghxx(indy,:);
        C(1:y_nbr , id_z6_xf_xf_xf) = 1/6*ghxxx(indy,:);

        D = zeros(y_nbr,e_nbr);
        D(1:y_nbr , id_e1_u      ) = ghu(indy,:) + 1/2*ghuss(indy,:);
        D(1:y_nbr , id_e2_u_u    ) = 1/2*ghuu(indy,:);
        D(1:y_nbr , id_e3_xf_u   ) = ghxu(indy,:);
        D(1:y_nbr , id_e5_xs_u   ) = ghxu(indy,:);
        D(1:y_nbr , id_e7_xf_xf_u) = 1/2*ghxxu(indy,:);
        D(1:y_nbr , id_e10_xf_u_u) = 1/2*ghxuu(indy,:);
        D(1:y_nbr , id_e13_u_u_u ) = 1/6*ghuuu(indy,:);

        c = zeros(z_nbr,1);
        c(id_z2_xs    , 1) = 1/2*ghs2(indx,:) + 1/2*ghuu(indx,:)*E_uu(:);
        c(id_z3_xf_xf , 1) = hu_hu*E_uu(:);

        d = zeros(y_nbr,1);
        d(1:y_nbr , 1) = 0.5*ghs2(indy,:) + 0.5*ghuu(indy,:)*E_uu(:);

        Varinov = zeros(e_nbr,e_nbr);
        Varinov(id_e1_u       , id_e1_u      ) = E_uu;
        Varinov(id_e1_u       , id_e5_xs_u   ) = kron(E_xs',E_uu);
        Varinov(id_e1_u       , id_e6_u_xs   ) = kron(E_uu,E_xs');
        Varinov(id_e1_u       , id_e7_xf_xf_u) = kron(E_xfxf(:)',E_uu);
        Varinov(id_e1_u       , id_e8_xf_u_xf) = kron(E_uu,E_xfxf(:)')*kron(K_xu,speye(x_nbr));
        Varinov(id_e1_u       , id_e9_u_xf_xf) = kron(E_uu,E_xfxf(:)');
        Varinov(id_e1_u       , id_e13_u_u_u ) = reshape(E_u_u_u_u,u_nbr,u_nbr^3);
        
        Varinov(id_e2_u_u     , id_e2_u_u    ) = E_uu_uu-E_uu(:)*transpose(E_uu(:));
        
        Varinov(id_e3_xf_u    , id_e3_xf_u   ) = E_xfxf_E_uu;
        Varinov(id_e3_xf_u    , id_e4_u_xf   ) = E_xfxf_E_uu*K_ux';
        Varinov(id_e3_xf_u    , id_e5_xs_u   ) = E_xfxs_E_uu;
        Varinov(id_e3_xf_u    , id_e6_u_xs   ) = E_xfxs_E_uu*K_ux';
        Varinov(id_e3_xf_u    , id_e7_xf_xf_u) = kron(reshape(E_xf_xf_xf,x_nbr,x_nbr^2), E_uu);
        Varinov(id_e3_xf_u    , id_e8_xf_u_xf) = kron(reshape(E_xf_xf_xf,x_nbr,x_nbr^2), E_uu)*kron(speye(x_nbr),K_ux)';
        Varinov(id_e3_xf_u    , id_e9_u_xf_xf) = K_xu*kron(E_uu, reshape(E_xf_xf_xf, x_nbr, x_nbr^2));
        
        Varinov(id_e4_u_xf    , id_e3_xf_u   ) = K_ux*E_xfxf_E_uu;
        Varinov(id_e4_u_xf    , id_e4_u_xf   ) = kron(E_uu,E_xfxf);
        Varinov(id_e4_u_xf    , id_e5_xs_u   ) = K_ux*kron(E_xfxs,E_uu);
        Varinov(id_e4_u_xf    , id_e6_u_xs   ) = kron(E_uu, E_xfxs);
        Varinov(id_e4_u_xf    , id_e7_xf_xf_u) = K_ux*kron(reshape(E_xf_xf_xf,x_nbr,x_nbr^2),E_uu);
        Varinov(id_e4_u_xf    , id_e8_xf_u_xf) = kron(E_uu,reshape(E_xf_xf_xf,x_nbr,x_nbr^2))*kron(K_xu,speye(x_nbr))';
        Varinov(id_e4_u_xf    , id_e9_u_xf_xf) = kron(E_uu,reshape(E_xf_xf_xf,x_nbr,x_nbr^2));
        
        Varinov(id_e5_xs_u    , id_e1_u      ) = kron(E_xs, E_uu);
        Varinov(id_e5_xs_u    , id_e3_xf_u   ) = kron(E_xsxf, E_uu);
        Varinov(id_e5_xs_u    , id_e4_u_xf   ) = kron(E_xsxf, E_uu)*K_ux';
        Varinov(id_e5_xs_u    , id_e5_xs_u   ) = kron(E_xsxs, E_uu);
        Varinov(id_e5_xs_u    , id_e6_u_xs   ) = kron(E_xsxs, E_uu)*K_ux';
        Varinov(id_e5_xs_u    , id_e7_xf_xf_u) = kron(reshape(E_xs_xf_xf,x_nbr,x_nbr^2),E_uu);
        Varinov(id_e5_xs_u    , id_e8_xf_u_xf) = kron(reshape(E_xs_xf_xf,x_nbr,x_nbr^2),E_uu)*kron(speye(x_nbr),K_ux)';
        Varinov(id_e5_xs_u    , id_e9_u_xf_xf) = K_xu*kron(E_uu,reshape(E_xs_xf_xf,x_nbr,x_nbr^2));
        Varinov(id_e5_xs_u    , id_e13_u_u_u ) = kron(E_xs,reshape(E_u_u_u_u,u_nbr,u_nbr^3));
        
        Varinov(id_e6_u_xs    , id_e1_u      ) = kron(E_uu,E_xs);
        Varinov(id_e6_u_xs    , id_e3_xf_u   ) = K_ux*kron(E_xsxf, E_uu);
        Varinov(id_e6_u_xs    , id_e4_u_xf   ) = kron(E_uu, E_xsxf);
        Varinov(id_e6_u_xs    , id_e5_xs_u   ) = K_ux*kron(E_xsxs,E_uu);
        Varinov(id_e6_u_xs    , id_e6_u_xs   ) = kron(E_uu, E_xsxs);
        Varinov(id_e6_u_xs    , id_e7_xf_xf_u) = K_ux*kron(reshape(E_xs_xf_xf,x_nbr,x_nbr^2), E_uu);
        Varinov(id_e6_u_xs    , id_e8_xf_u_xf) = kron(E_uu, reshape(E_xs_xf_xf,x_nbr,x_nbr^2))*kron(K_xu,speye(x_nbr))';
        Varinov(id_e6_u_xs    , id_e9_u_xf_xf) = kron(E_uu, reshape(E_xs_xf_xf,x_nbr,x_nbr^2));
        Varinov(id_e6_u_xs    , id_e13_u_u_u ) = K_ux*kron(E_xs,reshape(E_u_u_u_u,u_nbr,u_nbr^3));
        
        Varinov(id_e7_xf_xf_u , id_e1_u      ) = kron(E_xfxf(:),E_uu);
        Varinov(id_e7_xf_xf_u , id_e3_xf_u   ) = kron(reshape(E_xf_xf_xf,x_nbr^2,x_nbr),E_uu);
        Varinov(id_e7_xf_xf_u , id_e4_u_xf   ) = kron(reshape(E_xf_xf_xf,x_nbr^2,x_nbr),E_uu)*K_ux';
        Varinov(id_e7_xf_xf_u , id_e5_xs_u   ) = kron(reshape(E_xf_xf_xs,x_nbr^2,x_nbr),E_uu);
        Varinov(id_e7_xf_xf_u , id_e6_u_xs   ) = kron(reshape(E_xf_xf_xs,x_nbr^2,x_nbr),E_uu)*K_ux';
        Varinov(id_e7_xf_xf_u , id_e7_xf_xf_u) = kron(reshape(E_xf_xf_xf_xf,x_nbr^2,x_nbr^2),E_uu);
        Varinov(id_e7_xf_xf_u , id_e8_xf_u_xf) = kron(reshape(E_xf_xf_xf_xf,x_nbr^2,x_nbr^2),E_uu)*kron(speye(x_nbr),K_ux)';
        Varinov(id_e7_xf_xf_u , id_e9_u_xf_xf) = kron(speye(x_nbr),K_ux)*kron(K_ux,speye(x_nbr))*kron(E_uu, reshape(E_xf_xf_xf_xf,x_nbr^2,x_nbr^2));
        Varinov(id_e7_xf_xf_u , id_e13_u_u_u ) = kron(E_xfxf(:),reshape(E_u_u_u_u,u_nbr,u_nbr^3));
        
        Varinov(id_e8_xf_u_xf , id_e1_u      ) = kron(K_xu,speye(x_nbr))*kron(E_uu,E_xfxf(:));
        Varinov(id_e8_xf_u_xf , id_e3_xf_u   ) = kron(speye(x_nbr),K_xu)*kron(reshape(E_xf_xf_xf,x_nbr^2,x_nbr),E_uu);
        Varinov(id_e8_xf_u_xf , id_e4_u_xf   ) = kron(K_xu,speye(x_nbr))*kron(E_uu,reshape(E_xf_xf_xf,x_nbr^2,x_nbr));
        Varinov(id_e8_xf_u_xf , id_e5_xs_u   ) = kron(speye(x_nbr),K_ux)*kron(reshape(E_xf_xf_xs,x_nbr^2,x_nbr),E_uu);
        Varinov(id_e8_xf_u_xf , id_e6_u_xs   ) = kron(K_xu,speye(x_nbr))*kron(E_uu,reshape(E_xf_xf_xs,x_nbr^2,x_nbr));
        Varinov(id_e8_xf_u_xf , id_e7_xf_xf_u) = kron(speye(x_nbr),K_ux)*kron(reshape(E_xf_xf_xf_xf,x_nbr^2,x_nbr^2),E_uu);
        Varinov(id_e8_xf_u_xf , id_e8_xf_u_xf) = kron(K_xu,speye(x_nbr))*kron(E_uu,reshape(E_xf_xf_xf_xf,x_nbr^2,x_nbr^2))*kron(K_xu,speye(x_nbr))';
        Varinov(id_e8_xf_u_xf , id_e9_u_xf_xf) = kron(K_xu,speye(x_nbr))*kron(E_uu,reshape(E_xf_xf_xf_xf,x_nbr^2,x_nbr^2));
        Varinov(id_e8_xf_u_xf , id_e13_u_u_u ) = kron(K_xu,speye(x_nbr))*kron(reshape(E_u_u_u_u,u_nbr,u_nbr^3),E_xfxf(:));
        
        Varinov(id_e9_u_xf_xf , id_e1_u      ) = kron(E_uu, E_xfxf(:));
        Varinov(id_e9_u_xf_xf , id_e3_xf_u   ) = kron(E_uu, reshape(E_xf_xf_xf,x_nbr^2,x_nbr))*K_xu';
        Varinov(id_e9_u_xf_xf , id_e4_u_xf   ) = kron(E_uu, reshape(E_xf_xf_xf,x_nbr^2,x_nbr));
        Varinov(id_e9_u_xf_xf , id_e5_xs_u   ) = kron(E_uu, reshape(E_xf_xf_xs,x_nbr^2,x_nbr))*K_xu';
        Varinov(id_e9_u_xf_xf , id_e6_u_xs   ) = kron(E_uu, reshape(E_xf_xf_xs,x_nbr^2,x_nbr));
        Varinov(id_e9_u_xf_xf , id_e7_xf_xf_u) = kron(speye(x_nbr),K_ux)*kron(K_ux,speye(x_nbr))*kron(reshape(E_xf_xf_xf_xf,x_nbr^2,x_nbr^2),E_uu);
        Varinov(id_e9_u_xf_xf , id_e8_xf_u_xf) = kron(E_uu,reshape(E_xf_xf_xf_xf,x_nbr^2,x_nbr^2))*kron(speye(x_nbr),K_ux)';
        Varinov(id_e9_u_xf_xf , id_e9_u_xf_xf) = kron(E_uu,reshape(E_xf_xf_xf_xf,x_nbr^2,x_nbr^2));
        Varinov(id_e9_u_xf_xf , id_e13_u_u_u ) = kron(reshape(E_u_u_u_u,u_nbr,u_nbr^3),E_xfxf(:));
        
        Varinov(id_e10_xf_u_u , id_e10_xf_u_u) = kron(E_xfxf,reshape(E_u_u_u_u,u_nbr^2,u_nbr^2));
        Varinov(id_e10_xf_u_u , id_e11_u_xf_u) = kron(E_xfxf,reshape(E_u_u_u_u,u_nbr^2,u_nbr^2))*kron(K_ux,speye(u_nbr))';
        Varinov(id_e10_xf_u_u , id_e12_u_u_xf) = kron(E_xfxf,reshape(E_u_u_u_u,u_nbr^2,u_nbr^2))*kron(K_ux,speye(u_nbr))'*kron(speye(u_nbr),K_ux)';
        
        Varinov(id_e11_u_xf_u , id_e10_xf_u_u) = kron(K_ux,speye(u_nbr))*kron(E_xfxf,reshape(E_u_u_u_u,u_nbr^2,u_nbr^2));
        Varinov(id_e11_u_xf_u , id_e11_u_xf_u) = kron(K_ux,speye(u_nbr))*kron(E_xfxf,reshape(E_u_u_u_u,u_nbr^2,u_nbr^2))*kron(K_ux,speye(u_nbr))';
        Varinov(id_e11_u_xf_u , id_e12_u_u_xf) = kron(speye(u_nbr),K_ux)*kron(reshape(E_u_u_u_u,u_nbr^2,u_nbr^2),E_xfxf);
        
        Varinov(id_e12_u_u_xf , id_e10_xf_u_u) = kron(K_ux,speye(u_nbr))*kron(speye(u_nbr),K_ux)*kron(E_xfxf,reshape(E_u_u_u_u,u_nbr^2,u_nbr^2));
        Varinov(id_e12_u_u_xf , id_e11_u_xf_u) = kron(reshape(E_u_u_u_u,u_nbr^2,u_nbr^2),E_xfxf)*kron(speye(u_nbr),K_xu)';
        Varinov(id_e12_u_u_xf , id_e12_u_u_xf) = kron(reshape(E_u_u_u_u,u_nbr^2,u_nbr^2),E_xfxf);
        
        Varinov(id_e13_u_u_u , id_e1_u       ) = reshape(E_u_u_u_u,u_nbr^3,u_nbr);
        Varinov(id_e13_u_u_u , id_e5_xs_u    ) = kron(E_xs', reshape(E_u_u_u_u,u_nbr^3,u_nbr));
        Varinov(id_e13_u_u_u , id_e6_u_xs    ) = kron(reshape(E_u_u_u_u,u_nbr^3,u_nbr),E_xs');
        Varinov(id_e13_u_u_u , id_e7_xf_xf_u ) = kron(E_xfxf(:)',reshape(E_u_u_u_u,u_nbr^3,u_nbr));
        Varinov(id_e13_u_u_u , id_e8_xf_u_xf ) = kron(E_xfxf(:)',reshape(E_u_u_u_u,u_nbr^3,u_nbr))*kron(speye(x_nbr),K_ux)';
        Varinov(id_e13_u_u_u , id_e9_u_xf_xf ) = kron(reshape(E_u_u_u_u,u_nbr^3,u_nbr), E_xfxf(:)');
        Varinov(id_e13_u_u_u , id_e13_u_u_u  ) = reshape(Q6Pu*E_u_u_u_u_u_u,u_nbr^3,u_nbr^3);
        
        E_z = (speye(z_nbr)-A)\c;
        
        
        if compute_derivs
            dA = zeros(z_nbr,z_nbr,totparam_nbr);
            dB = zeros(z_nbr,e_nbr,totparam_nbr);
            dC = zeros(y_nbr,z_nbr,totparam_nbr);
            dD = zeros(y_nbr,e_nbr,totparam_nbr);
            dc = zeros(z_nbr,totparam_nbr);
            dd = zeros(y_nbr,totparam_nbr);
            dVarinov = zeros(e_nbr,e_nbr,totparam_nbr);
            dE_z =zeros(z_nbr,totparam_nbr);

            for jp=1:totparam_nbr
                dA(id_z1_xf       , id_z1_xf       ,jp) = dghx(indx,:,jp);
                dA(id_z2_xs       , id_z2_xs       ,jp) = dghx(indx,:,jp);
                dA(id_z2_xs       , id_z3_xf_xf    ,jp) = 1/2*dghxx(indx,:,jp);
                dA(id_z3_xf_xf    , id_z3_xf_xf    ,jp) = dhx_hx(:,:,jp);
                dA(id_z4_xrd      , id_z1_xf       ,jp) = 3/6*dghxss(indx,:,jp);
                dA(id_z4_xrd      , id_z4_xrd      ,jp) = dghx(indx,:,jp);
                dA(id_z4_xrd      , id_z5_xf_xs    ,jp) = dghxx(indx,:,jp);
                dA(id_z4_xrd      , id_z6_xf_xf_xf ,jp) = 1/6*dghxxx(indx,:,jp);
                dA(id_z5_xf_xs    , id_z1_xf       ,jp) = kron(dghx(indx,:,jp),1/2*ghs2(indx,:)) + kron(ghx(indx,:),1/2*dghs2(indx,jp));
                dA(id_z5_xf_xs    , id_z5_xf_xs    ,jp) = dhx_hx(:,:,jp);
                dA(id_z5_xf_xs    , id_z6_xf_xf_xf ,jp) = kron(dghx(indx,:,jp),1/2*ghxx(indx,:)) + kron(ghx(indx,:),1/2*dghxx(indx,:,jp));
                dA(id_z6_xf_xf_xf , id_z6_xf_xf_xf ,jp) = kron(dghx(indx,:,jp),hx_hx) + kron(ghx(indx,:),dhx_hx(:,:,jp));
                
                dB(id_z1_xf       , id_e1_u       , jp) = dghu(indx,:,jp);
                dB(id_z2_xs       , id_e2_u_u     , jp) = 1/2*dghuu(indx,:,jp);
                dB(id_z2_xs       , id_e3_xf_u    , jp) = dghxu(indx,:,jp);
                dB(id_z3_xf_xf    , id_e2_u_u     , jp) = dhu_hu(:,:,jp);
                dB(id_z3_xf_xf    , id_e3_xf_u    , jp) = dhx_hu(:,:,jp);
                dB(id_z3_xf_xf    , id_e4_u_xf    , jp) = dhu_hx(:,:,jp);
                dB(id_z4_xrd      , id_e1_u       , jp) = 3/6*dghuss(indx,:,jp);
                dB(id_z4_xrd      , id_e5_xs_u    , jp) = dghxu(indx,:,jp);
                dB(id_z4_xrd      , id_e7_xf_xf_u , jp) = 3/6*dghxxu(indx,:,jp);
                dB(id_z4_xrd      , id_e10_xf_u_u , jp) = 3/6*dghxuu(indx,:,jp);
                dB(id_z4_xrd      , id_e13_u_u_u  , jp) = 1/6*dghuuu(indx,:,jp);
                dB(id_z5_xf_xs    , id_e1_u       , jp) = kron(dghu(indx,:,jp),1/2*ghs2(indx,:)) + kron(ghu(indx,:),1/2*dghs2(indx,jp));
                dB(id_z5_xf_xs    , id_e6_u_xs    , jp) = dhu_hx(:,:,jp);
                dB(id_z5_xf_xs    , id_e7_xf_xf_u , jp) = kron(dghx(indx,:,jp),ghxu(indx,:)) + kron(ghx(indx,:),dghxu(indx,:,jp));
                dB(id_z5_xf_xs    , id_e9_u_xf_xf , jp) = kron(dghu(indx,:,jp),1/2*ghxx(indx,:)) + kron(ghu(indx,:),1/2*dghxx(indx,:,jp));
                dB(id_z5_xf_xs    , id_e10_xf_u_u , jp) = kron(dghx(indx,:,jp),1/2*ghuu(indx,:)) + kron(ghx(indx,:),1/2*dghuu(indx,:,jp));
                dB(id_z5_xf_xs    , id_e11_u_xf_u , jp) = kron(dghu(indx,:,jp),ghxu(indx,:)) + kron(ghu(indx,:),dghxu(indx,:,jp));
                dB(id_z5_xf_xs    , id_e13_u_u_u  , jp) = kron(dghu(indx,:,jp),1/2*ghuu(indx,:)) + kron(ghu(indx,:),1/2*dghuu(indx,:,jp));
                dB(id_z6_xf_xf_xf , id_e7_xf_xf_u , jp) = kron(dhx_hx(:,:,jp),ghu(indx,:)) + kron(hx_hx,dghu(indx,:,jp));
                dB(id_z6_xf_xf_xf , id_e8_xf_u_xf , jp) = kron(dghx(indx,:,jp),hu_hx) + kron(ghx(indx,:),dhu_hx(:,:,jp));
                dB(id_z6_xf_xf_xf , id_e9_u_xf_xf , jp) = kron(dghu(indx,:,jp),hx_hx) + kron(ghu(indx,:),dhx_hx(:,:,jp));
                dB(id_z6_xf_xf_xf , id_e10_xf_u_u , jp) = kron(dhx_hu(:,:,jp),ghu(indx,:)) + kron(hx_hu,dghu(indx,:,jp));
                dB(id_z6_xf_xf_xf , id_e11_u_xf_u , jp) = kron(dghu(indx,:,jp),hx_hu) + kron(ghu(indx,:),dhx_hu(:,:,jp));
                dB(id_z6_xf_xf_xf , id_e12_u_u_xf , jp) = kron(dhu_hu(:,:,jp),ghx(indx,:)) + kron(hu_hu,dghx(indx,:,jp));
                dB(id_z6_xf_xf_xf , id_e13_u_u_u  , jp) = kron(dghu(indx,:,jp),hu_hu) + kron(ghu(indx,:),dhu_hu(:,:,jp));
        
                dC(1:y_nbr , id_z1_xf       , jp) = dghx(indy,:,jp) + 1/2*dghxss(indy,:,jp);
                dC(1:y_nbr , id_z2_xs       , jp) = dghx(indy,:,jp);
                dC(1:y_nbr , id_z3_xf_xf    , jp) = 1/2*dghxx(indy,:,jp);
                dC(1:y_nbr , id_z4_xrd      , jp) = dghx(indy,:,jp);
                dC(1:y_nbr , id_z5_xf_xs    , jp) = dghxx(indy,:,jp);
                dC(1:y_nbr , id_z6_xf_xf_xf , jp) = 1/6*dghxxx(indy,:,jp);

                dD(1:y_nbr , id_e1_u       , jp) = dghu(indy,:,jp) + 1/2*dghuss(indy,:,jp);
                dD(1:y_nbr , id_e2_u_u     , jp) = 1/2*dghuu(indy,:,jp);
                dD(1:y_nbr , id_e3_xf_u    , jp) = dghxu(indy,:,jp);
                dD(1:y_nbr , id_e5_xs_u    , jp) = dghxu(indy,:,jp);
                dD(1:y_nbr , id_e7_xf_xf_u , jp) = 1/2*dghxxu(indy,:,jp);
                dD(1:y_nbr , id_e10_xf_u_u , jp) = 1/2*dghxuu(indy,:,jp);
                dD(1:y_nbr , id_e13_u_u_u  , jp) = 1/6*dghuuu(indy,:,jp);
        
                dc(id_z2_xs    , jp) = 1/2*dghs2(indx,jp) + 1/2*dghuu(indx,:,jp)*E_uu(:) + 1/2*ghuu(indx,:)*vec(dE_uu(:,:,jp));
                dc(id_z3_xf_xf , jp) = dhu_hu(:,:,jp)*E_uu(:) + hu_hu*vec(dE_uu(:,:,jp));

                dd(1:y_nbr , jp) = 0.5*dghs2(indy,jp) + 0.5*dghuu(indy,:,jp)*E_uu(:) + 0.5*ghuu(indy,:)*vec(dE_uu(:,:,jp));

                dVarinov(id_e1_u       , id_e1_u       , jp) = dE_uu(:,:,jp);
                dVarinov(id_e1_u       , id_e5_xs_u    , jp) = kron(dE_xs(:,jp)',E_uu) + kron(E_xs',dE_uu(:,:,jp));
                dVarinov(id_e1_u       , id_e6_u_xs    , jp) = kron(dE_uu(:,:,jp),E_xs') + kron(E_uu,dE_xs(:,jp)');
                dVarinov(id_e1_u       , id_e7_xf_xf_u , jp) = kron(vec(dE_xfxf(:,:,jp))',E_uu) + kron(E_xfxf(:)',dE_uu(:,:,jp));
                dVarinov(id_e1_u       , id_e8_xf_u_xf , jp) = (kron(dE_uu(:,:,jp),E_xfxf(:)') + kron(E_uu,vec(dE_xfxf(:,:,jp))'))*kron(K_xu,speye(x_nbr));
                dVarinov(id_e1_u       , id_e9_u_xf_xf , jp) = kron(dE_uu(:,:,jp),E_xfxf(:)') + kron(E_uu,vec(dE_xfxf(:,:,jp))');
                if jp <= (stderrparam_nbr+corrparam_nbr)
                    dVarinov(id_e1_u       , id_e13_u_u_u  , jp) = reshape(QPu*dE_u_u_u_u(:,jp),u_nbr,u_nbr^3);
                    dVarinov(id_e2_u_u     , id_e2_u_u     , jp) = reshape(QPu*dE_u_u_u_u(:,jp),u_nbr^2,u_nbr^2) - vec(dE_uu(:,:,jp))*transpose(E_uu(:)) - E_uu(:)*transpose(vec(dE_uu(:,:,jp)));
                end
                
                dVarinov(id_e3_xf_u    , id_e3_xf_u    , jp) = dE_xfxf_E_uu(:,:,jp);
                dVarinov(id_e3_xf_u    , id_e4_u_xf    , jp) = dE_xfxf_E_uu(:,:,jp)*K_ux';
                dVarinov(id_e3_xf_u    , id_e5_xs_u    , jp) = dE_xfxs_E_uu(:,:,jp);
                dVarinov(id_e3_xf_u    , id_e6_u_xs    , jp) = dE_xfxs_E_uu(:,:,jp)*K_ux';
                dVarinov(id_e3_xf_u    , id_e7_xf_xf_u , jp) = kron(reshape(dE_xf_xf_xf(:,jp),x_nbr,x_nbr^2), E_uu) + kron(reshape(E_xf_xf_xf,x_nbr,x_nbr^2), dE_uu(:,:,jp));
                dVarinov(id_e3_xf_u    , id_e8_xf_u_xf , jp) = (kron(reshape(dE_xf_xf_xf(:,jp),x_nbr,x_nbr^2), E_uu) + kron(reshape(E_xf_xf_xf,x_nbr,x_nbr^2), dE_uu(:,:,jp)))*kron(speye(x_nbr),K_ux)';
                dVarinov(id_e3_xf_u    , id_e9_u_xf_xf , jp) = K_xu*(kron(dE_uu(:,:,jp), reshape(E_xf_xf_xf, x_nbr, x_nbr^2)) + kron(E_uu, reshape(dE_xf_xf_xf(:,jp), x_nbr, x_nbr^2)));

                dVarinov(id_e4_u_xf    , id_e3_xf_u    , jp) = K_ux*dE_xfxf_E_uu(:,:,jp);
                dVarinov(id_e4_u_xf    , id_e4_u_xf    , jp) = kron(dE_uu(:,:,jp),E_xfxf) + kron(E_uu,dE_xfxf(:,:,jp));
                dVarinov(id_e4_u_xf    , id_e5_xs_u    , jp) = K_ux*(kron(dE_xfxs(:,:,jp),E_uu) + kron(E_xfxs,dE_uu(:,:,jp)));
                dVarinov(id_e4_u_xf    , id_e6_u_xs    , jp) = kron(dE_uu(:,:,jp), E_xfxs) + kron(E_uu, dE_xfxs(:,:,jp));
                dVarinov(id_e4_u_xf    , id_e7_xf_xf_u , jp) = K_ux*(kron(reshape(dE_xf_xf_xf(:,jp),x_nbr,x_nbr^2),E_uu) + kron(reshape(E_xf_xf_xf,x_nbr,x_nbr^2),dE_uu(:,:,jp)));
                dVarinov(id_e4_u_xf    , id_e8_xf_u_xf , jp) = (kron(dE_uu(:,:,jp),reshape(E_xf_xf_xf,x_nbr,x_nbr^2)) + kron(E_uu,reshape(dE_xf_xf_xf(:,jp),x_nbr,x_nbr^2)))*kron(K_xu,speye(x_nbr))';
                dVarinov(id_e4_u_xf    , id_e9_u_xf_xf , jp) = kron(dE_uu(:,:,jp),reshape(E_xf_xf_xf,x_nbr,x_nbr^2)) + kron(E_uu,reshape(dE_xf_xf_xf(:,jp),x_nbr,x_nbr^2));

                dVarinov(id_e5_xs_u    , id_e1_u       , jp) = kron(dE_xs(:,jp), E_uu) + kron(E_xs, dE_uu(:,:,jp));
                dVarinov(id_e5_xs_u    , id_e3_xf_u    , jp) = kron(dE_xsxf(:,:,jp), E_uu) + kron(E_xsxf, dE_uu(:,:,jp));
                dVarinov(id_e5_xs_u    , id_e4_u_xf    , jp) = (kron(dE_xsxf(:,:,jp), E_uu) + kron(E_xsxf, dE_uu(:,:,jp)))*K_ux';
                dVarinov(id_e5_xs_u    , id_e5_xs_u    , jp) = kron(dE_xsxs(:,:,jp), E_uu) + kron(E_xsxs, dE_uu(:,:,jp));
                dVarinov(id_e5_xs_u    , id_e6_u_xs    , jp) = (kron(dE_xsxs(:,:,jp), E_uu) + kron(E_xsxs, dE_uu(:,:,jp)))*K_ux';
                dVarinov(id_e5_xs_u    , id_e7_xf_xf_u , jp) = kron(reshape(dE_xs_xf_xf(:,jp),x_nbr,x_nbr^2),E_uu) + kron(reshape(E_xs_xf_xf,x_nbr,x_nbr^2),dE_uu(:,:,jp));
                dVarinov(id_e5_xs_u    , id_e8_xf_u_xf , jp) = (kron(reshape(dE_xs_xf_xf(:,jp),x_nbr,x_nbr^2),E_uu) + kron(reshape(E_xs_xf_xf,x_nbr,x_nbr^2),dE_uu(:,:,jp)))*kron(speye(x_nbr),K_ux)';
                dVarinov(id_e5_xs_u    , id_e9_u_xf_xf , jp) = K_xu*(kron(dE_uu(:,:,jp),reshape(E_xs_xf_xf,x_nbr,x_nbr^2)) + kron(E_uu,reshape(dE_xs_xf_xf(:,jp),x_nbr,x_nbr^2)));
                if jp <= (stderrparam_nbr+corrparam_nbr)
                    dVarinov(id_e5_xs_u    , id_e13_u_u_u  , jp) = kron(dE_xs(:,jp),reshape(E_u_u_u_u,u_nbr,u_nbr^3)) + kron(E_xs,reshape(QPu*dE_u_u_u_u(:,jp),u_nbr,u_nbr^3));
                else
                    dVarinov(id_e5_xs_u    , id_e13_u_u_u  , jp) = kron(dE_xs(:,jp),reshape(E_u_u_u_u,u_nbr,u_nbr^3));
                end

                dVarinov(id_e6_u_xs    , id_e1_u       , jp) = kron(dE_uu(:,:,jp),E_xs) + kron(E_uu,dE_xs(:,jp));
                dVarinov(id_e6_u_xs    , id_e3_xf_u    , jp) = K_ux*(kron(dE_xsxf(:,:,jp), E_uu) + kron(E_xsxf, dE_uu(:,:,jp)));
                dVarinov(id_e6_u_xs    , id_e4_u_xf    , jp) = kron(dE_uu(:,:,jp), E_xsxf) + kron(E_uu, dE_xsxf(:,:,jp));
                dVarinov(id_e6_u_xs    , id_e5_xs_u    , jp) = K_ux*(kron(dE_xsxs(:,:,jp),E_uu) + kron(E_xsxs,dE_uu(:,:,jp)));
                dVarinov(id_e6_u_xs    , id_e6_u_xs    , jp) = kron(dE_uu(:,:,jp), E_xsxs) + kron(E_uu, dE_xsxs(:,:,jp));
                dVarinov(id_e6_u_xs    , id_e7_xf_xf_u , jp) = K_ux*(kron(reshape(dE_xs_xf_xf(:,jp),x_nbr,x_nbr^2), E_uu) + kron(reshape(E_xs_xf_xf,x_nbr,x_nbr^2), dE_uu(:,:,jp)));
                dVarinov(id_e6_u_xs    , id_e8_xf_u_xf , jp) = (kron(dE_uu(:,:,jp), reshape(E_xs_xf_xf,x_nbr,x_nbr^2)) + kron(E_uu, reshape(dE_xs_xf_xf(:,jp),x_nbr,x_nbr^2)))*kron(K_xu,speye(x_nbr))';
                dVarinov(id_e6_u_xs    , id_e9_u_xf_xf , jp) = kron(dE_uu(:,:,jp), reshape(E_xs_xf_xf,x_nbr,x_nbr^2)) + kron(E_uu, reshape(dE_xs_xf_xf(:,jp),x_nbr,x_nbr^2));
                if jp <= (stderrparam_nbr+corrparam_nbr)
                    dVarinov(id_e6_u_xs    , id_e13_u_u_u  , jp) = K_ux*(kron(dE_xs(:,jp),reshape(E_u_u_u_u,u_nbr,u_nbr^3)) + kron(E_xs,reshape(QPu*dE_u_u_u_u(:,jp),u_nbr,u_nbr^3)));
                else
                    dVarinov(id_e6_u_xs    , id_e13_u_u_u  , jp) = K_ux*kron(dE_xs(:,jp),reshape(E_u_u_u_u,u_nbr,u_nbr^3));
                end

                dVarinov(id_e7_xf_xf_u , id_e1_u       , jp) = kron(vec(dE_xfxf(:,:,jp)),E_uu) + kron(E_xfxf(:),dE_uu(:,:,jp));
                dVarinov(id_e7_xf_xf_u , id_e3_xf_u    , jp) = kron(reshape(dE_xf_xf_xf(:,jp),x_nbr^2,x_nbr),E_uu) + kron(reshape(E_xf_xf_xf,x_nbr^2,x_nbr),dE_uu(:,:,jp));
                dVarinov(id_e7_xf_xf_u , id_e4_u_xf    , jp) = (kron(reshape(dE_xf_xf_xf(:,jp),x_nbr^2,x_nbr),E_uu) + kron(reshape(E_xf_xf_xf,x_nbr^2,x_nbr),dE_uu(:,:,jp)))*K_ux';
                dVarinov(id_e7_xf_xf_u , id_e5_xs_u    , jp) = kron(reshape(dE_xf_xf_xs(:,jp),x_nbr^2,x_nbr),E_uu) + kron(reshape(E_xf_xf_xs,x_nbr^2,x_nbr),dE_uu(:,:,jp));
                dVarinov(id_e7_xf_xf_u , id_e6_u_xs    , jp) = (kron(reshape(dE_xf_xf_xs(:,jp),x_nbr^2,x_nbr),E_uu) + kron(reshape(E_xf_xf_xs,x_nbr^2,x_nbr),dE_uu(:,:,jp)))*K_ux';
                dVarinov(id_e7_xf_xf_u , id_e7_xf_xf_u , jp) = kron(reshape(dE_xf_xf_xf_xf(:,jp),x_nbr^2,x_nbr^2),E_uu) + kron(reshape(E_xf_xf_xf_xf,x_nbr^2,x_nbr^2),dE_uu(:,:,jp));
                dVarinov(id_e7_xf_xf_u , id_e8_xf_u_xf , jp) = (kron(reshape(dE_xf_xf_xf_xf(:,jp),x_nbr^2,x_nbr^2),E_uu) + kron(reshape(E_xf_xf_xf_xf,x_nbr^2,x_nbr^2),dE_uu(:,:,jp)))*kron(speye(x_nbr),K_ux)';
                dVarinov(id_e7_xf_xf_u , id_e9_u_xf_xf , jp) = kron(speye(x_nbr),K_ux)*kron(K_ux,speye(x_nbr))*(kron(dE_uu(:,:,jp), reshape(E_xf_xf_xf_xf,x_nbr^2,x_nbr^2)) + kron(E_uu, reshape(dE_xf_xf_xf_xf(:,jp),x_nbr^2,x_nbr^2)));
                if jp <= (stderrparam_nbr+corrparam_nbr)
                    dVarinov(id_e7_xf_xf_u , id_e13_u_u_u  , jp) = kron(vec(dE_xfxf(:,:,jp)),reshape(E_u_u_u_u,u_nbr,u_nbr^3)) + kron(E_xfxf(:),reshape(QPu*dE_u_u_u_u(:,jp),u_nbr,u_nbr^3));
                else
                    dVarinov(id_e7_xf_xf_u , id_e13_u_u_u  , jp) = kron(vec(dE_xfxf(:,:,jp)),reshape(E_u_u_u_u,u_nbr,u_nbr^3));
                end

                dVarinov(id_e8_xf_u_xf , id_e1_u       , jp) = kron(K_xu,speye(x_nbr))*(kron(dE_uu(:,:,jp),E_xfxf(:)) + kron(E_uu,vec(dE_xfxf(:,:,jp))));
                dVarinov(id_e8_xf_u_xf , id_e3_xf_u    , jp) = kron(speye(x_nbr),K_xu)*(kron(reshape(dE_xf_xf_xf(:,jp),x_nbr^2,x_nbr),E_uu) + kron(reshape(E_xf_xf_xf,x_nbr^2,x_nbr),dE_uu(:,:,jp)));
                dVarinov(id_e8_xf_u_xf , id_e4_u_xf    , jp) = kron(K_xu,speye(x_nbr))*(kron(dE_uu(:,:,jp),reshape(E_xf_xf_xf,x_nbr^2,x_nbr)) + kron(E_uu,reshape(dE_xf_xf_xf(:,jp),x_nbr^2,x_nbr)));
                dVarinov(id_e8_xf_u_xf , id_e5_xs_u    , jp) = kron(speye(x_nbr),K_ux)*(kron(reshape(dE_xf_xf_xs(:,jp),x_nbr^2,x_nbr),E_uu) + kron(reshape(E_xf_xf_xs,x_nbr^2,x_nbr),dE_uu(:,:,jp)));
                dVarinov(id_e8_xf_u_xf , id_e6_u_xs    , jp) = kron(K_xu,speye(x_nbr))*(kron(dE_uu(:,:,jp),reshape(E_xf_xf_xs,x_nbr^2,x_nbr)) + kron(E_uu,reshape(dE_xf_xf_xs(:,jp),x_nbr^2,x_nbr)));
                dVarinov(id_e8_xf_u_xf , id_e7_xf_xf_u , jp) = kron(speye(x_nbr),K_ux)*(kron(reshape(dE_xf_xf_xf_xf(:,jp),x_nbr^2,x_nbr^2),E_uu) + kron(reshape(E_xf_xf_xf_xf,x_nbr^2,x_nbr^2),dE_uu(:,:,jp)));
                dVarinov(id_e8_xf_u_xf , id_e8_xf_u_xf , jp) = kron(K_xu,speye(x_nbr))*(kron(dE_uu(:,:,jp),reshape(E_xf_xf_xf_xf,x_nbr^2,x_nbr^2)) + kron(E_uu,reshape(dE_xf_xf_xf_xf(:,jp),x_nbr^2,x_nbr^2)))*kron(K_xu,speye(x_nbr))';
                dVarinov(id_e8_xf_u_xf , id_e9_u_xf_xf , jp) = kron(K_xu,speye(x_nbr))*(kron(dE_uu(:,:,jp),reshape(E_xf_xf_xf_xf,x_nbr^2,x_nbr^2)) + kron(E_uu,reshape(dE_xf_xf_xf_xf(:,jp),x_nbr^2,x_nbr^2)));
                if jp <= (stderrparam_nbr+corrparam_nbr)
                    dVarinov(id_e8_xf_u_xf , id_e13_u_u_u  , jp) = kron(K_xu,speye(x_nbr))*(kron(reshape(QPu*dE_u_u_u_u(:,jp),u_nbr,u_nbr^3),E_xfxf(:)) + kron(reshape(E_u_u_u_u,u_nbr,u_nbr^3),vec(dE_xfxf(:,:,jp))));
                else
                    dVarinov(id_e8_xf_u_xf , id_e13_u_u_u  , jp) = kron(K_xu,speye(x_nbr))*kron(reshape(E_u_u_u_u,u_nbr,u_nbr^3),vec(dE_xfxf(:,:,jp)));
                end

                dVarinov(id_e9_u_xf_xf , id_e1_u       , jp) = kron(dE_uu(:,:,jp), E_xfxf(:)) + kron(E_uu, vec(dE_xfxf(:,:,jp)));
                dVarinov(id_e9_u_xf_xf , id_e3_xf_u    , jp) = (kron(dE_uu(:,:,jp), reshape(E_xf_xf_xf,x_nbr^2,x_nbr)) + kron(E_uu, reshape(dE_xf_xf_xf(:,jp),x_nbr^2,x_nbr)))*K_xu';
                dVarinov(id_e9_u_xf_xf , id_e4_u_xf    , jp) = kron(dE_uu(:,:,jp), reshape(E_xf_xf_xf,x_nbr^2,x_nbr)) + kron(E_uu, reshape(dE_xf_xf_xf(:,jp),x_nbr^2,x_nbr));
                dVarinov(id_e9_u_xf_xf , id_e5_xs_u    , jp) = (kron(dE_uu(:,:,jp), reshape(E_xf_xf_xs,x_nbr^2,x_nbr)) + kron(E_uu, reshape(dE_xf_xf_xs(:,jp),x_nbr^2,x_nbr)))*K_xu';
                dVarinov(id_e9_u_xf_xf , id_e6_u_xs    , jp) = kron(dE_uu(:,:,jp), reshape(E_xf_xf_xs,x_nbr^2,x_nbr)) + kron(E_uu, reshape(dE_xf_xf_xs(:,jp),x_nbr^2,x_nbr));
                dVarinov(id_e9_u_xf_xf , id_e7_xf_xf_u , jp) = kron(speye(x_nbr),K_ux)*kron(K_ux,speye(x_nbr))*(kron(reshape(dE_xf_xf_xf_xf(:,jp),x_nbr^2,x_nbr^2),E_uu) + kron(reshape(E_xf_xf_xf_xf,x_nbr^2,x_nbr^2),dE_uu(:,:,jp)));
                dVarinov(id_e9_u_xf_xf , id_e8_xf_u_xf , jp) = (kron(dE_uu(:,:,jp),reshape(E_xf_xf_xf_xf,x_nbr^2,x_nbr^2)) + kron(E_uu,reshape(dE_xf_xf_xf_xf(:,jp),x_nbr^2,x_nbr^2)))*kron(speye(x_nbr),K_ux)';
                dVarinov(id_e9_u_xf_xf , id_e9_u_xf_xf , jp) = kron(dE_uu(:,:,jp),reshape(E_xf_xf_xf_xf,x_nbr^2,x_nbr^2)) + kron(E_uu,reshape(dE_xf_xf_xf_xf(:,jp),x_nbr^2,x_nbr^2));
                if jp <= (stderrparam_nbr+corrparam_nbr)
                    dVarinov(id_e9_u_xf_xf , id_e13_u_u_u  , jp) = kron(reshape(QPu*dE_u_u_u_u(:,jp),u_nbr,u_nbr^3),E_xfxf(:)) + kron(reshape(E_u_u_u_u,u_nbr,u_nbr^3),vec(dE_xfxf(:,:,jp)));
                else
                    dVarinov(id_e9_u_xf_xf , id_e13_u_u_u  , jp) = kron(reshape(E_u_u_u_u,u_nbr,u_nbr^3),vec(dE_xfxf(:,:,jp)));
                end
                
                if jp <= (stderrparam_nbr+corrparam_nbr)
                    dVarinov(id_e10_xf_u_u , id_e10_xf_u_u , jp) = kron(dE_xfxf(:,:,jp),reshape(E_u_u_u_u,u_nbr^2,u_nbr^2)) + kron(E_xfxf,reshape(QPu*dE_u_u_u_u(:,jp),u_nbr^2,u_nbr^2));
                    dVarinov(id_e10_xf_u_u , id_e11_u_xf_u , jp) = (kron(dE_xfxf(:,:,jp),reshape(E_u_u_u_u,u_nbr^2,u_nbr^2)) + kron(E_xfxf,reshape(QPu*dE_u_u_u_u(:,jp),u_nbr^2,u_nbr^2)))*kron(K_ux,speye(u_nbr))';
                    dVarinov(id_e10_xf_u_u , id_e12_u_u_xf , jp) = (kron(dE_xfxf(:,:,jp),reshape(E_u_u_u_u,u_nbr^2,u_nbr^2)) + kron(E_xfxf,reshape(QPu*dE_u_u_u_u(:,jp),u_nbr^2,u_nbr^2)))*kron(K_ux,speye(u_nbr))'*kron(speye(u_nbr),K_ux)';
                else
                    dVarinov(id_e10_xf_u_u , id_e10_xf_u_u , jp) = kron(dE_xfxf(:,:,jp),reshape(E_u_u_u_u,u_nbr^2,u_nbr^2));
                    dVarinov(id_e10_xf_u_u , id_e11_u_xf_u , jp) = kron(dE_xfxf(:,:,jp),reshape(E_u_u_u_u,u_nbr^2,u_nbr^2))*kron(K_ux,speye(u_nbr))';
                    dVarinov(id_e10_xf_u_u , id_e12_u_u_xf , jp) = kron(dE_xfxf(:,:,jp),reshape(E_u_u_u_u,u_nbr^2,u_nbr^2))*kron(K_ux,speye(u_nbr))'*kron(speye(u_nbr),K_ux)';
                end
                
                if jp <= (stderrparam_nbr+corrparam_nbr)
                    dVarinov(id_e11_u_xf_u , id_e10_xf_u_u , jp) = kron(K_ux,speye(u_nbr))*(kron(dE_xfxf(:,:,jp),reshape(E_u_u_u_u,u_nbr^2,u_nbr^2)) + kron(E_xfxf,reshape(QPu*dE_u_u_u_u(:,jp),u_nbr^2,u_nbr^2)));
                    dVarinov(id_e11_u_xf_u , id_e11_u_xf_u , jp) = kron(K_ux,speye(u_nbr))*(kron(dE_xfxf(:,:,jp),reshape(E_u_u_u_u,u_nbr^2,u_nbr^2)) + kron(E_xfxf,reshape(QPu*dE_u_u_u_u(:,jp),u_nbr^2,u_nbr^2)))*kron(K_ux,speye(u_nbr))';
                    dVarinov(id_e11_u_xf_u , id_e12_u_u_xf , jp) = kron(speye(u_nbr),K_ux)*(kron(reshape(QPu*dE_u_u_u_u(:,jp),u_nbr^2,u_nbr^2),E_xfxf) + kron(reshape(E_u_u_u_u,u_nbr^2,u_nbr^2),dE_xfxf(:,:,jp)));
                else
                    dVarinov(id_e11_u_xf_u , id_e10_xf_u_u , jp) = kron(K_ux,speye(u_nbr))*kron(dE_xfxf(:,:,jp),reshape(E_u_u_u_u,u_nbr^2,u_nbr^2));
                    dVarinov(id_e11_u_xf_u , id_e11_u_xf_u , jp) = kron(K_ux,speye(u_nbr))*kron(dE_xfxf(:,:,jp),reshape(E_u_u_u_u,u_nbr^2,u_nbr^2))*kron(K_ux,speye(u_nbr))';
                    dVarinov(id_e11_u_xf_u , id_e12_u_u_xf , jp) = kron(speye(u_nbr),K_ux)*kron(reshape(E_u_u_u_u,u_nbr^2,u_nbr^2),dE_xfxf(:,:,jp));
                end
                
                if jp <= (stderrparam_nbr+corrparam_nbr)
                    dVarinov(id_e12_u_u_xf , id_e10_xf_u_u , jp) = kron(K_ux,speye(u_nbr))*kron(speye(u_nbr),K_ux)*(kron(dE_xfxf(:,:,jp),reshape(E_u_u_u_u,u_nbr^2,u_nbr^2)) + kron(E_xfxf,reshape(QPu*dE_u_u_u_u(:,jp),u_nbr^2,u_nbr^2)));
                    dVarinov(id_e12_u_u_xf , id_e11_u_xf_u , jp) = (kron(reshape(QPu*dE_u_u_u_u(:,jp),u_nbr^2,u_nbr^2),E_xfxf) + kron(reshape(E_u_u_u_u,u_nbr^2,u_nbr^2),dE_xfxf(:,:,jp)))*kron(speye(u_nbr),K_xu)';
                    dVarinov(id_e12_u_u_xf , id_e12_u_u_xf , jp) = kron(reshape(QPu*dE_u_u_u_u(:,jp),u_nbr^2,u_nbr^2),E_xfxf) + kron(reshape(E_u_u_u_u,u_nbr^2,u_nbr^2),dE_xfxf(:,:,jp));
                else
                    dVarinov(id_e12_u_u_xf , id_e10_xf_u_u , jp) = kron(K_ux,speye(u_nbr))*kron(speye(u_nbr),K_ux)*kron(dE_xfxf(:,:,jp),reshape(E_u_u_u_u,u_nbr^2,u_nbr^2));
                    dVarinov(id_e12_u_u_xf , id_e11_u_xf_u , jp) = kron(reshape(E_u_u_u_u,u_nbr^2,u_nbr^2),dE_xfxf(:,:,jp))*kron(speye(u_nbr),K_xu)';
                    dVarinov(id_e12_u_u_xf , id_e12_u_u_xf , jp) = kron(reshape(E_u_u_u_u,u_nbr^2,u_nbr^2),dE_xfxf(:,:,jp));
                end
                
                if jp <= (stderrparam_nbr+corrparam_nbr)
                    dVarinov(id_e13_u_u_u , id_e1_u        , jp) = reshape(QPu*dE_u_u_u_u(:,jp),u_nbr^3,u_nbr);
                    dVarinov(id_e13_u_u_u , id_e5_xs_u     , jp) = kron(dE_xs(:,jp)', reshape(E_u_u_u_u,u_nbr^3,u_nbr)) + kron(E_xs', reshape(QPu*dE_u_u_u_u(:,jp),u_nbr^3,u_nbr));
                    dVarinov(id_e13_u_u_u , id_e6_u_xs     , jp) = kron(reshape(QPu*dE_u_u_u_u(:,jp),u_nbr^3,u_nbr),E_xs') + kron(reshape(E_u_u_u_u,u_nbr^3,u_nbr),dE_xs(:,jp)');
                    dVarinov(id_e13_u_u_u , id_e7_xf_xf_u  , jp) = kron(vec(dE_xfxf(:,:,jp))',reshape(E_u_u_u_u,u_nbr^3,u_nbr)) + kron(E_xfxf(:)',reshape(QPu*dE_u_u_u_u(:,jp),u_nbr^3,u_nbr));
                    dVarinov(id_e13_u_u_u , id_e8_xf_u_xf  , jp) = (kron(vec(dE_xfxf(:,:,jp))',reshape(E_u_u_u_u,u_nbr^3,u_nbr)) + kron(E_xfxf(:)',reshape(QPu*dE_u_u_u_u(:,jp),u_nbr^3,u_nbr)))*kron(speye(x_nbr),K_ux)';
                    dVarinov(id_e13_u_u_u , id_e9_u_xf_xf  , jp) = kron(reshape(QPu*dE_u_u_u_u(:,jp),u_nbr^3,u_nbr), E_xfxf(:)') + kron(reshape(E_u_u_u_u,u_nbr^3,u_nbr), vec(dE_xfxf(:,:,jp))');
                else                    
                    dVarinov(id_e13_u_u_u , id_e5_xs_u     , jp) = kron(dE_xs(:,jp)', reshape(E_u_u_u_u,u_nbr^3,u_nbr));
                    dVarinov(id_e13_u_u_u , id_e6_u_xs     , jp) = kron(reshape(E_u_u_u_u,u_nbr^3,u_nbr),dE_xs(:,jp)');
                    dVarinov(id_e13_u_u_u , id_e7_xf_xf_u  , jp) = kron(vec(dE_xfxf(:,:,jp))',reshape(E_u_u_u_u,u_nbr^3,u_nbr));
                    dVarinov(id_e13_u_u_u , id_e8_xf_u_xf  , jp) = kron(vec(dE_xfxf(:,:,jp))',reshape(E_u_u_u_u,u_nbr^3,u_nbr))*kron(speye(x_nbr),K_ux)';
                    dVarinov(id_e13_u_u_u , id_e9_u_xf_xf  , jp) = kron(reshape(E_u_u_u_u,u_nbr^3,u_nbr), vec(dE_xfxf(:,:,jp))');
                end
                if jp <= (stderrparam_nbr+corrparam_nbr)
                    dVarinov(id_e13_u_u_u , id_e13_u_u_u   , jp) = reshape(Q6Pu*dE_u_u_u_u_u_u(:,jp),u_nbr^3,u_nbr^3);
                end

                dE_z(:,jp) = (speye(z_nbr)-A)\(dc(:,jp) + dA(:,:,jp)*E_z);
            end
        end
    end
end


E_y = Yss(indy,:) + C*E_z + d;
Om_z = B*Varinov*transpose(B);
Om_y = D*Varinov*transpose(D);

if compute_derivs
    dE_y = zeros(y_nbr,totparam_nbr);
    dOm_z = zeros(z_nbr,z_nbr,totparam_nbr);
    dOm_y = zeros(y_nbr,y_nbr,totparam_nbr);
    for jp = 1:totparam_nbr
        dE_y(:,jp)   = dC(:,:,jp)*E_z + C*dE_z(:,jp) + dd(:,jp);
        dOm_z(:,:,jp) = dB(:,:,jp)*Varinov*B' + B*dVarinov(:,:,jp)*B' + B*Varinov*dB(:,:,jp)';
        dOm_y(:,:,jp) = dD(:,:,jp)*Varinov*D' + D*dVarinov(:,:,jp)*D' + D*Varinov*dD(:,:,jp)';
        if jp > (stderrparam_nbr+corrparam_nbr)
            dE_y(:,jp) = dE_y(:,jp) + dYss(indy,jp-stderrparam_nbr-corrparam_nbr); %add steady state
        end
    end
end


%% Store into output structure
dr.pruned.indx = indx;
dr.pruned.indy = indy;
%dr.pruned.E_xfxf = E_xfxf;
dr.pruned.A = A;
dr.pruned.B = B;
dr.pruned.C = C;
dr.pruned.D = D;
dr.pruned.c = c;
dr.pruned.d = d;
dr.pruned.Om_z = Om_z;
dr.pruned.Om_y = Om_y;
dr.pruned.Varinov = Varinov;
dr.pruned.E_z = E_z;
dr.pruned.E_y = E_y;
if compute_derivs == 1
    %dr.pruned.dE_xfxf = dE_xfxf;
    dr.pruned.dA = dA;
    dr.pruned.dB = dB;
    dr.pruned.dC = dC;
    dr.pruned.dD = dD;
    dr.pruned.dc = dc;
    dr.pruned.dd = dd;
    dr.pruned.dOm_z = dOm_z;
    dr.pruned.dOm_y = dOm_y;
    dr.pruned.dVarinov = dVarinov;
    dr.pruned.dE_z = dE_z;
    dr.pruned.dE_y = dE_y;
end
