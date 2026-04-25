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
    disp('Errox\re nelle date!')
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

        % state at time t
        X_t = State(1:2,ind1);
        Z= Z(1:17,1:2,[ind1,ind2]);
        A=d(1:17,[ind1,ind2]);
        
        a = [1;0];
        DZ = squeeze(Z(:,:,2)) - squeeze(Z(:,:,1));
        DA = A(:,2) - A(:,1);
        Z = squeeze(Z(:,:,2));

    case 'AFNS'

        lambda = result.params.lambda;
        K = result.params.KP;
        barX = result.params.thetaP;
        Sigma_X = result.params.Sigma;
        X_t = result.filteredStates(ind1,:)';
        Z_fun = @(tau)[ones(size(tau)), (1-exp(-lambda*tau))./(lambda*tau), (1-exp(-lambda*tau))./(lambda*tau) - exp(-lambda*tau)];
 
        a = [1;1;0];
        DZ = zeros(length(tau_mat),length(X_t));
        DA = 0;
        Z = Z_fun(tau_mat');

end 


[EY,alpha,beta]=EYgivenDeltaF2(a,tau_J,h,X_t,K,barX,Sigma_X, Z, DZ, DA, Delta_F);

out.Ey=EY;
out.alpha=alpha;
out.beta=beta;


end



function [EDY_cond,EDY,beta] = EYgivenDeltaF2( ...
    a,...                         % specify the linear comnìbinations of factors that gives the short rate
    T1,h, ...                     % tempi
    X_t, ...                      % stato al tempo t
    K,barX, Sigma_X, ...          % parametri OU
    Z, DZ, DA, ...                % funzioni di carico
    DeltaF )                      % valore osservato ΔF_t
    
    % EYgivenDeltaF2  :  E_t[ Δ_h y(t,T) | Δ_h F_t ]
    
    S = Sigma_X*Sigma_X';
    Lambda_h = expm(-K*h);
    mu_h = (eye(length(a)) - Lambda_h)*barX;
    
    Lambda_fun = @(u)(expm(-K*u));
    Q_h = integral(@(u) Lambda_fun(u)*S*Lambda_fun(u)', ...
                              0, h, 'ArrayValued', true, 'RelTol',1e-10,'AbsTol',1e-12);
    
    b = Lambda_fun(T1-h)' * a;

    beta = (Z*Q_h*b) / (b.'*Q_h*b);      % slope
    EDY = DA + DZ*X_t + Z*((Lambda_h -eye(length(a)))*X_t + mu_h);
    EDY_cond   = EDY + beta * DeltaF;                  % Et[Δy | ΔF]

end