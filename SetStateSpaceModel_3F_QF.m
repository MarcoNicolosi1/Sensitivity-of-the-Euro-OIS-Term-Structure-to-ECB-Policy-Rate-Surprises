function [Z, d, GG, T_mtx, c, HH, a, P] = SetStateSpaceModel_3F_QF(hyper_par, t, T, t_J, N_Estr)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[n_obs, N] = size(T); %maturities, N = number of series for estimation
m_states = 3;     % number of states
g = N+m_states;   % number of errors


% Measurement equation
% y_{t} = Z_t \alpha_{t} +  d_t + G_t \epsilon_t

% Transition equation
% \alpha_{t+1} = T_mtx_t \alpha_{t} + c_t+ H_t \epsilon_t

Z = nan(N_Estr,m_states,n_obs);
d = nan(N_Estr,n_obs);

Z = nan(N,m_states,n_obs);
d = nan(N,n_obs);

%a = nan(n,N);
%b = nan(n,N);
T_mtx = nan(m_states,m_states,n_obs);
c = nan(m_states,n_obs);
HH = nan(m_states,m_states,n_obs);



% gJ=parJ{2};
% wJ=parJ{3};
% 

k_xi = exp(hyper_par(1));
sig_xi = exp(hyper_par(2));
k_theta = exp(hyper_par(3));
theta_bar = hyper_par(4);
sig_theta = exp(hyper_par(5));




par1{1} = k_xi;
par1{2} = sig_xi;
par2{1} = k_theta;
par2{2} = theta_bar;
par2{3} = sig_theta;
parJ{1} = 0;%aJ
parJ{2} = [1; 0];% gJ;
parJ{3} = exp(hyper_par(6));% w_J % std dev of jumps
w = parJ{3}; 
rho = 0.9*(1 - 2*exp(hyper_par(7))/(1+exp(hyper_par(7))));% correlazione tra le due variabili di stato xi e theta

sigma = exp(hyper_par(8));  % std dev measurement error 

if length(hyper_par) == 8
    lambda_xi = 0; % market price of risk of xi
    lambda_theta = 0; % market price of risk of theta
    gamma_J = 0;
else
    lambda_xi = hyper_par(9); % market price of risk of xi
    lambda_theta = hyper_par(10); % market price of risk of theta
    gamma_J = exp(hyper_par(11)); % market price of risk for jumps
end

div = 360;

for i = 1:n_obs
    for j = 1:N_Estr
      %  target_t = ESTR(i);
        target_t = 0;        
        [a, b] = YieldBond_Affine_2f(t(i),T(i,j), t_J, par1,par2,rho,parJ, target_t);
        Z(j,:,i) = [b' 1];
        d(j,i) = a;
    end

    if i == n_obs
        dt = days(t(n_obs)-t(n_obs-1))/div;
    else
        dt = days(t(i+1)-t(i))/div;
    end
    
    X_bar_Obj = -[-1/k_xi -1/k_theta; 0 -1/k_theta]*[sig_xi*lambda_xi; k_theta*theta_bar + sig_theta*(rho*lambda_xi + sqrt(1-rho^2)*lambda_theta)];
   
    A_Obj = A_Obj_2f(dt,k_theta,k_xi,lambda_xi,lambda_theta,rho,sig_xi,sig_theta,theta_bar);
   
    B_Obj = B_Obj_2f(dt,k_theta,k_xi);
    V_Obj = V_Obj_2f(dt,k_theta,k_xi,rho,sig_xi,sig_theta);

    d_J = sum(t(i) == t_J);

    T_mtx(:,:,i) = [B_Obj, [0; 0]; [gamma_J*d_J 0 1]];    
        
    c(:,i)  = [A_Obj; 0];

    cholV = chol(V_Obj,"lower");
    HH(:,:,i) = [zeros(m_states,N), [cholV, [0;0]; [0 0 w*d_J]] ]*[zeros(m_states,N), [cholV, [0;0]; [0 0 w*d_J]]]'; % H_t * H_t'  

end

GG = [sigma*eye(N,N) zeros(N,m_states)]*[sigma*eye(N,N) zeros(N,m_states)]'; % G_t * G_t' 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Initial conditions 
a = [X_bar_Obj; 0];     % a_{1|0}

tau = 1000;
P = [V_Obj_2f(tau,k_theta,k_xi,rho,sig_xi,sig_theta) [0; 0]; [0 0 w^2]]; % P_{1|0}     

end


