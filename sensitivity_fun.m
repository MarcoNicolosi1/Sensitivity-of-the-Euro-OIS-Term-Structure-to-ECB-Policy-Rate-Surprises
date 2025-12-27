function out=sensitivity_fun(t_0,Delta_F,h_shift,file_upload)

    % Compute E_t0 [dy(t0,T_vec)|Delta_F]
    
    load(file_upload)
    kxi  = k_xi;   kth  = k_theta;
    barT = theta_bar;
       sxi  = sig_xi;  sth = sig_theta;
    gamma= gamma_J;
    
    ind1=find(t_y==t_0); % 
    ind2=ind1+h_shift;
    
    h  = years(t_y(ind2)-t_0);
    if isempty(h) 
        disp('Error on dates!')
        EY=[];
    else    
        % jump time (first scheduled)
        T1 = min(t_J(t_J>t_0));
        tau_J=years(T1-t_0); % time to jump in years
        %indJ=find(t_y==t_J);
        
        
        % state at time t
        xi_t    = State(1,ind1);
        theta_t = State(2,ind2);
        
        
           
        Zxi=squeeze(Z(1:17,1,[ind1,ind2]));
        Zth=squeeze(Z(1:17,2,[ind1,ind2]));
        d_estr=d(1:17,[ind1,ind2]);
        
        Tvec=[1/365,2/365,7/365,1/12,3/12,6/12,9/12,1:10];
        Zxi_fun = @(s) interp1(Tvec,Zxi,s,'linear'); 
        Zth_fun = @(s) interp1(Tvec,Zth,s,'linear');
        d_fun   = @(s) interp1(Tvec,d_estr,s,'linear');
        
        
        tau_mat=Tvec;
        
        
        a1=gamma;
        a2=0;
        [EY,alpha,beta] = EYgivenDeltaF( ...
                      tau_mat,tau_J,h, ...              % times
                      xi_t,theta_t, ...                 % states at time t
                      kxi,kth,barT, ...                 % OU parameters
                      rho,sxi,sth, ...                  % correlation and volatilites
                      a1,a2,w_J, ...                    % jumps: a=[a1;a2], var_eps = w^2
                      d_fun,Zxi_fun,Zth_fun, ...        % factor loadings
                      Delta_F);
    end
    out.Ey=EY;
    out.alpha=alpha;
    out.beta=beta;
end



function [EY,c_det,beta] = EYgivenDeltaF( ...
              T,T1,h, ...                       % times
              xi_t,theta_t, ...                 % states at time t
              kxi,kth,barT, ...                 % OU parameters
              rho,sxi,sth, ...                  % correlation and volatilities
              a1,a2,w, ...                      % jumps: a=[a1;a2], var eps = w^2
              d_fun,Zxi_fun,Zth_fun, ...        % factor loadings
              DeltaF )                          % observed value ΔF_t

    % EYgivenDeltaF  :  E_t[ Δ_h y(t,T) | Δ_h F_t ]
    %
    % implement:
    %  - if T1>t+h:
    %      Δy = c + Z'ΔX ,   ΔF = b'ΔX,  β = (Z'Σ_h b)/(b'Σ_h b)
    %      EY = c + β*ΔF
    %  - if T1≤t+h:
    %      Δy = c + Z'ΔX + J , ΔF = J - F_t
    %      EY = c + F_t + (1 + Z'Cov(ΔX,J)/Var(J)) * ΔF
    %
    % where Σ_h = Var_t(ΔX),  b = Φ(τ-h)' a,  τ = T1-t,
    % Cov(ΔX,J) = [∫_0^τ Φ(h-u) S Φ(τ-u)' du] a,  Var(J)= a'V_T1 a + w^2,
    % V_T1 = ∫_0^τ Φ(s) S Φ(s)' ds,  a=[a1;a2].
    
    t=0; 
    
    a   = [a1; a2];
    S   = [ sxi^2, rho*sxi*sth; rho*sxi*sth, sth^2 ];
    eta = [0; kth*barT];
    
    tmp=Zxi_fun(T);
    Zxi_tT=tmp(:,1);
    Zxi_thT=tmp(:,2);
    
    tmp=Zth_fun(T);
    Zth_tT= tmp(:,1);
    Zth_thT=tmp(:,2);
    
    Z  = [Zxi_thT'; Zth_thT'];
    
    % deterministic part (known in t)
    tmp=d_fun(T);
    Delta_d     = tmp(:,2) - tmp(:,1);
    DeltaZshift = (Zxi_thT - Zxi_tT)*xi_t + (Zth_thT - Zth_tT)*theta_t;
    c_det       = Delta_d + DeltaZshift;
    
    tau = T1-t;   % distance between t=0 and a jump
    
    if T1 > t+h
        % ======= CASE A: no jumps in (t, t+h] =======
        % b = Φ(τ-h)' a
        b = subPhi(tau-h,kxi,kth) .' * a;
    
        % Σ_h = ∫_0^h Φ(s) S Φ(s)' ds
        Sigma_h = integral(@(s) subPhi(s,kxi,kth)*S*subPhi(s,kxi,kth).', ...
                           0, h, 'ArrayValued', true, 'RelTol',1e-10,'AbsTol',1e-12);
    
        beta = (Z.'*Sigma_h*b) / (b.'*Sigma_h*b);      % slope
        EY   = c_det + beta * DeltaF;                  % Et[Δy | ΔF]
    
    else
        % ======= CASE B: jump within (t, t+h] =======
        % Cov(ΔX, X_{T1}) = ∫_0^τ Φ(h-u) S Φ(τ-u)' du
        Cov_dX_XT1 = integral(@(u) subPhi(h-u,kxi,kth)*S*subPhi(tau-u,kxi,kth).', ...
                              0, tau, 'ArrayValued', true, 'RelTol',1e-10,'AbsTol',1e-12);
    
        % Cov(ΔX, J) = Cov(ΔX, X_{T1}) * a
        cov_dX_J = Cov_dX_XT1 * a;                     % 2x1
    
        % Var(J) = a' Var(X_{T1}) a + w^2,  Var(X_{T1}) = ∫_0^τ Φ(s) S Φ(s)' ds
        V_XT1 = integral(@(s) subPhi(s,kxi,kth)*S*subPhi(s,kxi,kth).', ...
                         0, tau, 'ArrayValued', true, 'RelTol',1e-10,'AbsTol',1e-12);
        varJ  = a.'*V_XT1*a + w^2;
    
        % F_t = a' μ_t ,   μ_t = Φ(τ) X_t + ∫_0^τ Φ(s) η ds
        mu_t  = subPhi(tau,kxi,kth)*[xi_t; theta_t] + ...
                integral(@(s) subPhi(s,kxi,kth)*eta, 0, tau, ...
                         'ArrayValued', true, 'RelTol',1e-10,'AbsTol',1e-12);
        Ft    = a.'*mu_t;
    
        beta  = 1 + (Z.'*cov_dX_J)/varJ;
        %EY    = c_det + Ft + beta * DeltaF;
        c_det = c_det + Ft;
        EY   = c_det + beta * DeltaF;                  % Et[Δy | ΔF]
    end
end

% ---------- local functions ----------
function Phi = subPhi(s,kx,kt)
% transition e^{K s} per K = [-kx, kx; 0, -kt]
if abs(kt-kx) > 1e-12
    Phi = [ exp(-kx*s),              kx/(kt-kx)*(exp(-kx*s)-exp(-kt*s));
            0         ,              exp(-kt*s) ];
else
    % limit kt -> kx  : off-diagonal = kx*s*exp(-kx*s)
    Phi = [ exp(-kx*s),  kx*s*exp(-kx*s);
            0         ,  exp(-kx*s)    ];
end
end