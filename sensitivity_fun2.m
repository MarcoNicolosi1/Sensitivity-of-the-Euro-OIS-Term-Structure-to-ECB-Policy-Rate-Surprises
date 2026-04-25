function out=sensitivity_fun2(t_0,Delta_F,h_shift,file_upload,model)
% calcola E_t0 [dy(t0,T_vec)|Delta_F]

%load dataset
load(file_upload)

if nargin == 4
    model = '3F';
end

ind1=find(t_y==t_0); % 
ind2=ind1+h_shift;

h  = years(t_y(ind2)-t_0);

if isempty(h) 
    display('Errox\re nelle date!')
    EY=[];
    return;
end

% jump time (first scheduled)
T1 = min(t_J(t_J>t_0));
tau_J=years(T1-t_0); % time to jump in years
%indJ=find(t_y==t_J);
Tvec=[1/365,2/365,7/365,1/12,3/12,6/12,9/12,1:10];

tau_mat=Tvec;

switch model
    case '3F'
        K = [k_xi -k_xi; 0 k_theta];
        barX = [theta_bar;theta_bar];
        Sigma_X  = [sig_xi 0;  rho*sig_theta sqrt(1-rho^2)*sig_theta];
        gamma= gamma_J;
        % state at time t
        X_t = State(1:2,ind1);
        Z= Z(1:17,1:2,[ind1,ind2]);
        A=d(1:17,[ind1,ind2]);
        
        [EY,alpha,beta]=EYgivenDeltaF_DJTA(tau_J,h,X_t,K,barX,Sigma_X,A ,Z, Delta_F);

    case 'AFNS'
        lambda = result.params.lambda;
        K = result.params.KP;
        barX = result.params.thetaP;
        Sigma_X = result.params.Sigma;
        X_t = result.filteredStates(ind1,:)';
        Z_fun = @(tau)[ones(size(tau)), (1-exp(-lambda*tau))./(lambda*tau), (1-exp(-lambda*tau))./(lambda*tau) - exp(-lambda*tau)];
  
        [EY,alpha,beta]=EYgivenDeltaF_AFNS(tau_mat,tau_J,h,X_t,K,barX,Sigma_X,Z_fun,Delta_F);

end 

out.Ey=EY;
out.alpha=alpha;
out.beta=beta;


end

function [EDY_cond,EDY,beta] = EYgivenDeltaF_AFNS( ...
              T,T1,h, ...                     % tempi
              X_t, ...                 % stato al tempo t
              K,barX, Sigma_X, ...                 % parametri OU          
              Z_fun, ...        % funzioni di carico
              DeltaF )                          % valore osservato ΔF_t
    % EYgivenDeltaF2  :  E_t[ Δ_h y(t,T) | Δ_h F_t ]
    %
    % Implementa:
    %  
    %      Δy =  Z'ΔX ,   ΔF = b'ΔX,  β = (Z'Σ_h b)/(b'Σ_h b)
    %      EY = c + β*ΔF
    %  - se T1≤t+h:
    %      Δy = c + Z'ΔX + J , ΔF = J - F_t
    %      EY = c + F_t + (1 + Z'Cov(ΔX,J)/Var(J)) * ΔF
    %
    % dove Σ_h = Var_t(ΔX),  b = Φ(τ-h)' a,  τ = T1-t,
    % Cov(ΔX,J) = [∫_0^τ Φ(h-u) S Φ(τ-u)' du] a,  Var(J)= a'V_T1 a + w^2,
    % V_T1 = ∫_0^τ Φ(s) S Φ(s)' ds,  a=[a1;a2].
    
    % ----- setup costanti e funzioni locali
    t=0; % ho messo qui l'origine!
    
    a = [1;1;0];
    
    S = Sigma_X*Sigma_X';
    Lambda_h = expm(-K*h);
    mu_h = (eye(length(a)) - Lambda_h)*barX;
    
    Lambda_fun = @(u)(expm(-K*u));
    Q_h = integral(@(u) Lambda_fun(u)*S*Lambda_fun(u)', ...
                              0, h, 'ArrayValued', true, 'RelTol',1e-10,'AbsTol',1e-12);
    
    b = Lambda_fun(T1-h)' * a;
    Z = Z_fun(T');
    beta = (Z*Q_h*b) / (b.'*Q_h*b);      % slope
    EDY = Z*((Lambda_h -eye(length(a)))*X_t + mu_h);
    EDY_cond   = EDY + beta * DeltaF;                  % Et[Δy | ΔF]

end

function [EDY_cond,EDY,beta] = EYgivenDeltaF_DJTA( ...
              T1,h, ...                     % tempi
              X_t, ...                 % stato al tempo t
              K,barX, Sigma_X, ...                 % parametri OU
              A,...
              Z, ...        % funzioni di carico
              DeltaF )                          % valore osservato ΔF_t
    % EYgivenDeltaF2  :  E_t[ Δ_h y(t,T) | Δ_h F_t ]
    %
    % Implementa:
    %  
    %      Δy =  Z'ΔX ,   ΔF = b'ΔX,  β = (Z'Σ_h b)/(b'Σ_h b)
    %      EY = c + β*ΔF
    %  - se T1≤t+h:
    %      Δy = c + Z'ΔX + J , ΔF = J - F_t
    %      EY = c + F_t + (1 + Z'Cov(ΔX,J)/Var(J)) * ΔF
    %
    % dove Σ_h = Var_t(ΔX),  b = Φ(τ-h)' a,  τ = T1-t,
    % Cov(ΔX,J) = [∫_0^τ Φ(h-u) S Φ(τ-u)' du] a,  Var(J)= a'V_T1 a + w^2,
    % V_T1 = ∫_0^τ Φ(s) S Φ(s)' ds,  a=[a1;a2].
    
    % ----- setup costanti e funzioni locali
    t=0; % ho messo qui l'origine!
    
    a = [1;0];
    
    S = Sigma_X*Sigma_X';
    Lambda_h = expm(-K*h);
    mu_h = (eye(length(a)) - Lambda_h)*barX;
    
    Lambda_fun = @(u)(expm(-K*u));
    Q_h = integral(@(u) Lambda_fun(u)*S*Lambda_fun(u)', ...
                              0, h, 'ArrayValued', true, 'RelTol',1e-10,'AbsTol',1e-12);
    
    b = Lambda_fun(T1-h)' * a;

    DZ = squeeze(Z(:,:,2)) - squeeze(Z(:,:,1));
    DA = A(:,2) - A(:,1);
    beta = (squeeze(Z(:,:,2))*Q_h*b) / (b.'*Q_h*b);      % slope
    EDY = DA + DZ*X_t + squeeze(Z(:,:,2))*((Lambda_h -eye(length(a)))*X_t + mu_h);
    EDY_cond   = EDY + beta * DeltaF;                  % Et[Δy | ΔF]

end