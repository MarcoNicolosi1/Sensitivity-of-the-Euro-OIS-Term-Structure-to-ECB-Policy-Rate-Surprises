clear 
close all
load data.mat

%% set the dataaset with the calibrated model. 
%% This is the output of scrpt_calibration
calibrated_model = 'calibrated_model.mat';



% This script identifies the dates of large moves in the OIS rate over the next maintenance period (MP).
% It then computes the corresponding zero-yields responses and the associated model-implied responses. 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% read the time series of the OIS over the next maintenance period
Delta_OIS.M1=10^2*diff(OISZCMP1.mid);
Delta_OIS.date=OISZCMP1.Date(2:end);

ind=isnan(Delta_OIS.M1);
Delta_OIS.M1(ind)=[];
Delta_OIS.date(ind)=[];

% Almost all changes above 0.1 occur exactly four days after an ECB meeting.
% For this reason, differences corresponding to ECB meeting − 3 and + 4 days
% are set to NaN to exclude non-meaningful movements.

time_interval=[-1,6];
Delta_OIS=clean(Delta_OIS,BCE_meetings,time_interval);


%% Identify the dates when the OIS over the next maintenance period experiences large moves
pct_d=5;
Y_d = prctile(Delta_OIS.M1,pct_d);
DJ_d.F1=Delta_OIS.M1(Delta_OIS.M1<Y_d); % negative surprises (top pct%)
DJ_d.date=Delta_OIS.date(Delta_OIS.M1<Y_d);
N_d=length(DJ_d.F1);


pct_u=95;
Y_u = prctile(Delta_OIS.M1,pct_u);
DJ_u.F1=Delta_OIS.M1(Delta_OIS.M1>Y_u); % positive surprises (top pct%)
DJ_u.date=Delta_OIS.date(Delta_OIS.M1>Y_u);
N_u=length(DJ_u.F1);


%% load the time series of Ester-zero-yields 
load(calibrated_model,'yEster')
load(calibrated_model,'t_y')

Delta_Ester.M1=10^4*diff(yEster);% in bps
Delta_Ester.date=t_y(2:end);
Delta_Ester=clean(Delta_Ester,BCE_meetings,time_interval);
Delta_Ester.y=Delta_Ester.M1;
Delta_Ester.M1=[];


%% Compute the observed variations and the model variations of the Ester-zero-yields on the dates with larger surprises
h_shift=1;

[Deltay_d,hat_Deltay_d]=compute_Delta(Delta_Ester,DJ_d,h_shift,calibrated_model);
[Deltay_u,hat_Deltay_u]=compute_Delta(Delta_Ester,DJ_u,h_shift,calibrated_model);


%% TABLE 1
PCT_y=prctile([Deltay_d.y;Deltay_u.y],[5,25,50,75,95]);
PCT_MP=prctile([DJ_d.F1;DJ_u.F1],[5,25,50,75,95]);
TABLE_PCT=[PCT_MP',PCT_y(:,[4,8,end])]';
disp('Table 1')
disp(TABLE_PCT)


%% Perform the regression in Euqation (1) 
ind=1;
for mat=1:17
    X=[DJ_u.F1];
    X1=[ones(size(X)),X];
    [b1u(:,ind),b1u_int(:,:,ind),b1u_r(:,ind) ] = regress([Deltay_u.y(:,mat)],X1);
    X=[DJ_d.F1];
    X1=[ones(size(X)),X];
    [b1d(:,ind),b1d_int(:,:,ind),b1d_r(:,ind) ] = regress([Deltay_d.y(:,mat)],X1);
    X=[DJ_d.F1;DJ_u.F1];
    X1=[ones(size(X)),X];
    [b1(:,ind),b1_int(:,:,ind),b1_r(:,ind) ] = regress([Deltay_d.y(:,mat);Deltay_u.y(:,mat)],X1);
    ind=ind+1;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% FIGURES

maturities={'1d','2d','1w','1m','3m','6m','9m','1y','2y','3y','4y','5y','6y','7y','8y','9y','10y'};

% FIGURE 1(a) of the paper

[DJ,ind]=sort(DJ_u.F1,'descend');
dd=DJ_u.date(ind);

figure();
for i=1:2
    t_0=dd(i);
    out1=sensitivity_fun(t_0,10^-4*DJ(i),h_shift,calibrated_model);
    subplot(1,2,i)
    if not(isempty(out1.Ey))
        plot([1:17],out1.Ey*10^4,'+-',[1:17],Delta_Ester.y(Delta_Ester.date==t_0,:),'o:' ,'linewidth',2)
        legend('model','realized','Location','southeast' )
    end
    title(['t = ', datestr(t_0),'; \Delta f_t =  ', num2str(DJ(i)), ' bps'])
    xticks([1:17])
    xticklabels({'1d','2d','1w','1m','3m','6m','9m','1y','2y','3y','4y','5y','6y','7y','8y','9y','10y'});
    xlim([1,17])
    ylim([-25 25])
end


%FIGURE 1(b) of the paper

[DJ,ind]=sort(DJ_d.F1,'ascend');
dd=DJ_d.date(ind);

figure();
for i=1:2
    t_0=dd(i);
    out1=sensitivity_fun(t_0,10^-4*DJ(i),h_shift,calibrated_model);
    subplot(1,2,i)
    plot([1:17],10^4*out1.Ey,'+-',[1:17],Delta_Ester.y(Delta_Ester.date==t_0,:),'o:' ,'linewidth',2)
    legend('model','realized','location',' northeast' )
    title(['t = ', datestr(t_0),'; \Delta f_t =  ', num2str(DJ(i)), ' bps'])
    xticks([1:17])
    xticklabels({'1d','2d','1w','1m','3m','6m','9m','1y','2y','3y','4y','5y','6y','7y','8y','9y','10y'});
    xlim([1,17])
    ylim([-25 25])
end


% FIGURE 2 of the paper
figure();
mat_b1=length(b1(2,:));
plot((1:mat_b1),b1(2,:),'o:','LineWidth',2)
xticks([1:mat_b1])
xticklabels(maturities(end-mat_b1+1:end));
xlim([1,17])
hold on
Beta_mod=[hat_Deltay_d.beta;hat_Deltay_d.beta];
M_beta_mod=mean(Beta_mod,1);
plot((1:17),M_beta_mod,'+:','LineWidth',2)
plot((1:mat_b1),squeeze(b1_int(2,:,:)),'.:','LineWidth',2)
grid on
legend('\beta(\tau)','mean(\beta^M(t,\tau))','  conf. int. \beta(\tau)','conf. int. \beta(\tau)')


% FIGURE 4 of the paper
figure();
start=1;
PCT=prctile([hat_Deltay_d.res;hat_Deltay_u.res],[5,25,50,75,95]);
plot((start:17),PCT(:,start:end),'+:','LineWidth',2)
xticks([start:17])
xticklabels(maturities(start:end));
xlim([start,17])
grid on
hold off


return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function Delta_OIS=clean(Delta_OIS,BCE_meetings,time_interval)
    for i=1:length(BCE_meetings)
        ind=and(Delta_OIS.date>=BCE_meetings(i)+days(time_interval(1)),Delta_OIS.date<=BCE_meetings(i)+days(time_interval(2)));
        Delta_OIS.M1(ind,:)=[];
        Delta_OIS.date(ind)=[];
    end
end

function  [Deltay,Deltay_hat]=compute_Delta(Delta_Ester,DJ,h_shift,calibrated_model)
    N=length(DJ.date);
    Deltay.y=NaN(N,17);%
    Deltay_hat.alpha=NaN(N,17);
    Deltay_hat.beta=NaN(N,17);
    Deltay_hat.y=NaN(N,17);
    for i=1:N
        ind=Delta_Ester.date==DJ.date(i);
        if sum(ind)>0
            Deltay.y(i,:)=Delta_Ester.y(ind,:);
            Deltay.date(i)=Delta_Ester.date(ind);
            Deltay_hat.date(i)=Delta_Ester.date(ind);
            out1=sensitivity_fun(DJ.date(i),10^-4*DJ.F1(i),h_shift,calibrated_model);
            Deltay_hat.y(i,:)=10^4*out1.Ey; %Ey=alpha+beta*DJ.F1
            Deltay_hat.alpha(i,:)=out1.alpha;
            Deltay_hat.beta(i,:)=out1.beta;
        end
        Deltay_hat.res=Deltay.y-Deltay_hat.y; % in bP
    end
end

function [RMSE_d_t,SS_d_t,sign_d_t,Rel_RMSE]=stats(Deltay,hat_Deltay)
    ERR_d=(Deltay-hat_Deltay); 
    Ncol=size(ERR_d,2);
    ERR_d=ERR_d(not(isnan(ERR_d)));
    ERR_d=reshape(ERR_d,[],Ncol);
    rel_ERR_d=ERR_d./Deltay;
    RMSE_d_t=sqrt(mean(ERR_d.^2,1));
    Rel_RMSE=sqrt(mean(rel_ERR_d.^2,1));
    
    
    Prod=Deltay.*hat_Deltay;
    Prod=Prod(not(isnan(Prod)));
    Prod=reshape(Prod,[],Ncol);
    DD=hat_Deltay(not(isnan(hat_Deltay)));
    DD=reshape(DD,[],Ncol);
    SS_d_t=sum(Prod,1)./sum(DD.^2,1);
    
    S=sign(Deltay);
    S_hat=sign(hat_Deltay);
    Concord=S.*S_hat;
    sign_d_t=sum(Concord==1,1)/size(Concord,1);
end
