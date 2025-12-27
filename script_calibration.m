clear all 
close all

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Script to estimate a 3-dimensional latent factor model that 
% captures the dynamics of the Estr OIS term structure.
%
% Let y_Estr(t,T) be the yield to maturity T observed at time t derived by the
% term structure of Estr OIS.
%
% y_Estr(t,T) = -1/(T-t)ln(P_Estr(t,T)) where
% P_Estr(t,T) = E_t exp(-\int_t^T r_s ds)
%
% r_t is a short latent rate that is constant up to the next BCE meeting
% date when it can jump. The jump size is normal with expectation \xi_t and
% variance w_J^2. 
% \xi_t is a latent factor that follows an Ornstein Oulenbeck process whose 
% long term drift \theta is another latent factor, following an Ornstein Oulenbeck
% process.
%
% r_t = r_0 + \sum_{i \in times of jumps} J_i
% J_t \sim N(\xi_t,w_J^2)
% d\xi_t = k_x*(\theta_t-\xi_t)*dt + sig_xi*dW1_t
% d\theta_t = k_x*(bar_theta-\theta_t)*dt + sig_theta*dW2_t
%
% The model is affine. Hence
% y_Estr = a + b*[\xi_t;\theta_t;r_t] 
%

%
% The model is estimated via maximum Log likelihood and Kalman filter
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

load data.mat

% Convert tables to timetable for data synchronization
OISZC = table2timetable(OISZC);

%Remove missing values
EstrZC = rmmissing(OISZC);

% Extract and normalize interest rates
EstrZY=[EstrZC.EURESTROISZCONPOINTZEROYIELD, EstrZC.EURESTROISZCTNPOINTZEROYIELD, EstrZC.EURESTROISZC1WPOINTZEROYIELD,...
    EstrZC.EURESTROISZC1MPOINTZEROYIELD,EstrZC.EURESTROISZC3MPOINTZEROYIELD, EstrZC.EURESTROISZC6MPOINTZEROYIELD, EstrZC.EURESTROISZC9MPOINTZEROYIELD,...    
    EstrZC.EURESTROISZC1YPOINTZEROYIELD, EstrZC.EURESTROISZC2YPOINTZEROYIELD, EstrZC.EURESTROISZC3YPOINTZEROYIELD, EstrZC.EURESTROISZC4YPOINTZEROYIELD, ...
    EstrZC.EURESTROISZC5YPOINTZEROYIELD,EstrZC.EURESTROISZC6YPOINTZEROYIELD, EstrZC.EURESTROISZC7YPOINTZEROYIELD, EstrZC.EURESTROISZC8YPOINTZEROYIELD,...
    EstrZC.EURESTROISZC9YPOINTZEROYIELD,EstrZC.EURESTROISZC10YPOINTZEROYIELD]/100;

[~, N_mat] = size(EstrZY);

% Configure maturities (same maturities for each obervation time)
t = EstrZC.Date;
matZC=repmat(t,1,N_mat);
matZC(:,1:2)=t+days(1:2);
matZC(:,3)=t+days(7);
m_shift=[1,3,6,9];
matZC(:,4)=datetime(year(t),month(t)+m_shift(1),day(t));
matZC(:,5)=datetime(year(t),month(t)+m_shift(2),day(t));
matZC(:,6)=datetime(year(t),month(t)+m_shift(3),day(t));
matZC(:,7)=datetime(year(t),month(t)+m_shift(4),day(t));
y_shift=1:10;
for i = 1:10
    matZC(:,8+(i-1)) = datetime(year(t)+y_shift(i),month(t),day(t));
end

div=360;
mat=days(matZC-t)./div;% maturities in years act/360
   

%% choose the time windows for estimation
t0 = 1; 
t_end = length(t);


%All observed maturities 
mat_Estr = mat(t0:t_end, 1:N_mat);
matZC_Estr = matZC(t0:t_end, 1:N_mat);
yEster =  log(1+EstrZY(t0:t_end, 1:N_mat));

%% Choose Maturities on Estr OIS used to calibrate the model
mat_in_Estr_cal = 1;
mat_end_Estr_cal = N_mat;
mat_Estr_cal = mat(t0:t_end, mat_in_Estr_cal:mat_end_Estr_cal);
matZC_Estr_cal = matZC(t0:t_end,mat_in_Estr_cal:mat_end_Estr_cal);
yEster_cal =  log(1+EstrZY(t0:t_end,mat_in_Estr_cal:mat_end_Estr_cal));

t_y = EstrZC.Date(t0:t_end);

%Dates of BCE meetings. Jumps may occur at BCE meetings plus 6 days
t_J = BCE_meetings+days(6);
t_J = [t_J; t_J(end)+ [42:42:42*80]']; 

[~, N_mat_Estr_cal] = size(yEster_cal);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Estimation

% Hyperparameters:   
% k_xi = exp(hyper_par(1));
% sig_xi = exp(hyper_par(2));
% k_theta = exp(hyper_par(3));
% theta_bar = hyper_par(4);
% sig_theta = exp(hyper_par(5));
% w_J = exp(hyper_par(6));% w_J % std dev of jumps
% rho = 1 - 2*exp(hyper_par(7))/(1+exp(hyper_par(7)));% correlation between xi and theta
% sigma = exp(hyper_par(8));  % std dev measurement error 
% lambda_xi = hyper_par(9); % market price of risk of xi
% lambda_theta = hyper_par(10); % market price of risk of theta
% gamma = exp(hyper_par(11)); % change of measure fo jumps
%


ProbMeas = 'Change of meausure';
%ProbMeas = 'RN';

% Change of measure
switch ProbMeas
    case 'RN'
        hyper_par0 = [0.4258; -4.7781; -0.1263; 0.0009;  -5.7466; -4.9906;  10.8253; -7.7348]; 
    case 'Change of meausure'  
        hyper_par0 = [0.4083; -5.4323; 0.3059; 0.0002;  -5.7018; -6.1787;  10.8253; -7.5671; 0; 0; 0];
end

%% Estimate parameters using maximum likelihood
f  = @(hyper_par)LogLik_fn(hyper_par, yEster_cal, t_y, matZC_Estr_cal, t_J, N_mat_Estr_cal); 
opts = optimset('Display','iter','TolX',1e-9,'TolFun',1e-9,...
                'Diagnostics','on', 'MaxIter',1000, 'MaxFunEvals', 10000,...
                'LargeScale', 'off', 'PlotFcns', @optimplotfval);
[hyper_par, fval, exitflag, output, grad, hessian] = fminunc(f, hyper_par0, opts);
 

disp(['MLE of transformed parameters:  ', num2str(hyper_par')]);
disp(['Likelihood:  ', num2str(-1.0*fval)]);
k_xi = exp(hyper_par(1));
sig_xi = exp(hyper_par(2));
k_theta = exp(hyper_par(3));
theta_bar = hyper_par(4);
sig_theta = exp(hyper_par(5));
w_J = exp(hyper_par(6));% w_J % std dev of jumps
rho = 0.9*(1 - 2*exp(hyper_par(7))/(1+exp(hyper_par(7))));
sigma = exp(hyper_par(8));  % std dev measurement error 


cov_par = inv(hessian);
sig_CI = sqrt(abs(diag(cov_par)));
CI_k_xi = k_xi*[1-1.96*sig_CI(1) 1+1.96*sig_CI(1)];
CI_sig_xi = sig_xi*[1-1.96*sig_CI(2) 1+1.96*sig_CI(2)];
CI_k_theta = k_theta*[1-1.96*sig_CI(3) 1+1.96*sig_CI(3)];
CI_theta_bar = [theta_bar-1.96*sig_CI(4) theta_bar+1.96*sig_CI(4)];
CI_sig_theta = sig_theta*[1-1.96*sig_CI(5) 1+1.96*sig_CI(5)];
CI_w_J = w_J*[1-1.96*sig_CI(6) 1+1.96*sig_CI(6)];
tmp = abs(-2*exp(hyper_par(7))/(1+exp(hyper_par(7)))^2);
CI_rho = [rho-1.96*tmp*sig_CI(7) rho+1.96*tmp*sig_CI(7)];
CI_sigma = sigma*[1-1.96*sig_CI(8) 1+1.96*sig_CI(8)];


if length(hyper_par) == 8
    lambda_xi = 0; % market price of risk of xi
    lambda_theta = 0; % market price of risk of theta
    gamma_J = 1;
else
    lambda_xi = hyper_par(9); % market price of risk of xi
    lambda_theta = hyper_par(10); % market price of risk of theta
    gamma_J = exp(hyper_par(11)); % market price of risk for jumps
    CI_lambda_xi = [lambda_xi-1.96*sig_CI(9) lambda_xi+1.96*sig_CI(9)];
    CI_lambda_theta = [lambda_theta-1.96*sig_CI(10) lambda_theta+1.96*sig_CI(10)];
    CI_gamma_J = gamma_J*[1-1.96*sig_CI(11) 1+1.96*sig_CI(11)];
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% display parameter estimates of TABLE 2

disp('ESTIMATION RESULTS         ');
disp(['k_xi = ', num2str(k_xi), '     CI: ', num2str(CI_k_xi)]);
disp(['sigma_xi = ', num2str(sig_xi), '     CI: ', num2str(CI_sig_xi)]);
disp(['k_theta = ', num2str(k_theta), '     CI: ', num2str(CI_k_theta)]);
disp(['theta_bar = ', num2str(theta_bar), '     CI: ', num2str(CI_theta_bar)]);
disp(['sigma_theta = ', num2str(sig_theta), '     CI: ', num2str(CI_sig_theta)]);
disp(['w_J = ', num2str(w_J), '     CI: ', num2str(CI_w_J)]);
disp(['rho = ', num2str(rho), '     CI: ', num2str(CI_rho)]);
disp(['sigma = ', num2str(sigma), '     CI: ', num2str(CI_sigma)]);

if length(hyper_par) == 11
    disp(['lambda_xi ', num2str(lambda_xi), '     CI: ', num2str(CI_lambda_xi)]);
    disp(['lambda_theta ', num2str(lambda_theta), '     CI: ', num2str(CI_lambda_theta)]);
    disp(['gamma_J ', num2str(gamma_J), '     CI: ', num2str(CI_gamma_J)]);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Create and estimate the state-space model
[Z, d, GG, T, c, HH, a, P] = SetStateSpaceModel_3F_QF(hyper_par, t_y, matZC_Estr, t_J, N_mat);

[Inn, F, StatePred, CovStatePred, LogLik] = ...
    KalmanFilter_fn(yEster(:,mat_in_Estr_cal:end), Z(mat_in_Estr_cal:end,:,:), ...
    d(mat_in_Estr_cal:end,:), GG(mat_in_Estr_cal:end,mat_in_Estr_cal:end), T, c, HH, a, P);

%% Smoothed estimates of latent states
[State, VarState] = Smoother_fn(yEster(:,mat_in_Estr_cal:end), Z(mat_in_Estr_cal:end,:,:), d(mat_in_Estr_cal:end,:), ...
    GG(mat_in_Estr_cal:end,mat_in_Estr_cal:end), T, c, HH, a, P);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% model term structure

yEster_model = yEster;
Nobs = length(t_y);

for i = 1:Nobs
    yEster_model(i,1:N_mat) = d(1:N_mat,i) + Z(1:N_mat,:,i)*State(:,i); 
end

%% RMSE
SE_Ester = (yEster_model-yEster).^2;
RSE_Ester = sqrt(mean(SE_Ester,2));
RMSE_Ester = sqrt(mean(SE_Ester))*10000;


% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% save results
% save('calibrated_model.mat')
% 


% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Figure

% FIGURE 3 of the paper
figure();
subplot(2,2,1)
i = 4;
plot(t_y,yEster_model(:,i))
hold on;
plot(t_y,yEster(:,i),'.','MarkerSize',3);
hold off;
legend('model','data','location','NorthWest')
title([num2str(mat_Estr_cal(1,i)*div/30,'%.0f'),' month yield'])
subplot(2,2,2)
i = 8;
plot(t_y,yEster_model(:,i))
hold on;
plot(t_y,yEster(:,i),'.','MarkerSize',3);
hold off;
title([num2str(mat_Estr_cal(1,i),'%.0f'),' year yield'])
subplot(2,2,3)
i = 10;
plot(t_y,yEster_model(:,i))
hold on;
plot(t_y,yEster(:,i),'.','MarkerSize',3);
hold off;
title([num2str(mat_Estr(1,i),'%.0f'),' year yield'])
subplot(2,2,4)
i = 17;
plot(t_y,yEster_model(:,i))
hold on;
plot(t_y,yEster(:,i),'.','MarkerSize',3);
hold off;
title([num2str(mat_Estr(1,i),'%.0f'),' year yield'])


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Other fugures
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% fit at a particular date
t_val = length(t_y);
figure();
plot(mat_Estr(t_val,:), 100*yEster(t_val,:),'o');
hold on
plot(mat_Estr(t_val,:), 100*yEster_model(t_val,:))
hold off
legend('market','model')
xlabel('maturities (years)')
title(['Model fit at ', datestr(t_y(t_val))])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% time-series fit for differente maturities 
figure();
sgtitle('Estimate for Estr')
subplot(3,3,1)
i = 1;
plot(t_y,yEster_model(:,i))
hold on;
plot(t_y,yEster(:,i),'.r');
hold off;
title(['yield at maturity = ', num2str(mat_Estr_cal(1,i)*div,'%.0f'),' days'])
subplot(3,3,2)
i = 3;
plot(t_y,yEster_model(:,i))
hold on;
plot(t_y,yEster(:,i),'.r');
hold off;
title(['yield at maturity = ', num2str(mat_Estr_cal(1,i)*div,'%.0f'),' days'])
subplot(3,3,3)
i = 5;
plot(t_y,yEster_model(:,i))
hold on;
plot(t_y,yEster(:,i),'.r');
hold off;
title(['yield at maturity = ', num2str(mat_Estr(1,i),'%.2f'),' years'])
subplot(3,3,4)
i = 6;
plot(t_y,yEster_model(:,i))
hold on;
plot(t_y,yEster(:,i),'.r');
hold off;
title(['yield at maturity = ', num2str(mat_Estr(1,i),'%.1f'),' years'])
subplot(3,3,5)
i = 8;
plot(t_y,yEster_model(:,i))
hold on;
plot(t_y,yEster(:,i),'.r');
hold off;
title(['yield at maturity = ', num2str(mat_Estr(1,i),'%.0f'),' year'])
subplot(3,3,6)
i = 9;
plot(t_y,yEster_model(:,i))
hold on;
plot(t_y,yEster(:,i),'.r');
hold off;
title(['yield at maturity = ', num2str(mat_Estr(1,i),'%.0f'),' years'])
subplot(3,3,7)
i = 12;
plot(t_y,yEster_model(:,i))
hold on;
plot(t_y,yEster(:,i),'.r');
hold off;
title(['yield at maturity = ', num2str(mat_Estr(1,i),'%.0f'),' years'])
subplot(3,3,8)
i = 14;
plot(t_y,yEster_model(:,i))
hold on;
plot(t_y,yEster(:,i),'.r');
hold off;
title(['yield at maturity = ', num2str(mat_Estr(1,i),'%.0f'),' years'])
subplot(3,3,9)
i = 17;
plot(t_y,yEster_model(:,i))
hold on;
plot(t_y,yEster(:,i),'.r');
hold off;
title(['yield at maturity = ', num2str(mat_Estr(1,i),'%.0f'),' years'])




