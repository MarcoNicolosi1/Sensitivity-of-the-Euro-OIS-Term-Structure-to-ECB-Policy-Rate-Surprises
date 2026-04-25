% AFNS calibration via Kalman-filter maximum likelihood
% Christensen, Diebold, and Rudebusch (2011)
%
% The implementation follows the paper's continuous-time AFNS setup under
% the P-measure, discretizes it exactly with matrix exponentials, and
% evaluates the Gaussian likelihood with the Kalman filter.
%
% The code supports:
%   - `modelType = 'correlated'`: unrestricted K^P and lower-triangular Sigma
%   - `modelType = 'independent'`: diagonal K^P and diagonal Sigma


clearvars;
clc;

load 'data_QF.mat'

% Convert tables to timetable for data synchronization
OISZC = table2timetable(OISZC);
OISZCMP1 = table2timetable(OISZCMP1);

% Synchronize and remove missing values
OISZC = synchronize(OISZCMP1, OISZC);          
EstrZC = rmmissing(OISZC);

% Extract and normalize interest rates
EstrZY=[EstrZC.EURESTROISZCONPOINTZEROYIELD, EstrZC.EURESTROISZCTNPOINTZEROYIELD, EstrZC.EURESTROISZC1WPOINTZEROYIELD,...
    EstrZC.EURESTROISZC1MPOINTZEROYIELD,EstrZC.EURESTROISZC3MPOINTZEROYIELD, EstrZC.EURESTROISZC6MPOINTZEROYIELD, EstrZC.EURESTROISZC9MPOINTZEROYIELD,...    
    EstrZC.EURESTROISZC1YPOINTZEROYIELD, EstrZC.EURESTROISZC2YPOINTZEROYIELD, EstrZC.EURESTROISZC3YPOINTZEROYIELD, EstrZC.EURESTROISZC4YPOINTZEROYIELD, ...
    EstrZC.EURESTROISZC5YPOINTZEROYIELD,EstrZC.EURESTROISZC6YPOINTZEROYIELD, EstrZC.EURESTROISZC7YPOINTZEROYIELD, EstrZC.EURESTROISZC8YPOINTZEROYIELD,...
    EstrZC.EURESTROISZC9YPOINTZEROYIELD,EstrZC.EURESTROISZC10YPOINTZEROYIELD]/100;
OIS_MP1 = OISZC.mid/100;
[~, N_mat] = size(EstrZY);




% Configure maturities (same maturities for each obervation time)
lag=0; %temporal lag
t = EstrZC.Date;
start=t+lag;
matZC=repmat(start,1,N_mat);
matZC(:,1:2)=start+days(1:2);
matZC(:,3)=start+days(7);
m_shift=[1,3,6,9];
matZC(:,4)=datetime(year(start),month(start)+m_shift(1),day(start));
matZC(:,5)=datetime(year(start),month(start)+m_shift(2),day(start));
matZC(:,6)=datetime(year(start),month(start)+m_shift(3),day(start));
matZC(:,7)=datetime(year(start),month(start)+m_shift(4),day(start));
y_shift=1:10;
for i = 1:10
    matZC(:,8+(i-1)) = datetime(year(start)+y_shift(i),month(start),day(start));
end

div=360;
mat=days(matZC-t)./div;% maturities in years act/360
   

%choose the time windows for estimation

t0 = 1; %all period
t_end = length(t);


%All observed maturities 
mat_Estr = mat(t0:t_end, 1:N_mat);
matZC_Estr = matZC(t0:t_end, 1:N_mat);
yEster =  log(1+EstrZY(t0:t_end, 1:N_mat));

%maturities on Estr OIS used to calibrate the model
mat_in_Estr_cal = 1;
mat_end_Estr_cal = N_mat;
mat_Estr_cal = mat(t0:t_end, mat_in_Estr_cal:mat_end_Estr_cal);
matZC_Estr_cal = matZC(t0:t_end,mat_in_Estr_cal:mat_end_Estr_cal);
yEster_cal =  log(1+EstrZY(t0:t_end,mat_in_Estr_cal:mat_end_Estr_cal));

t_y = EstrZC.Date(t0:t_end);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ----------------------------- User inputs ------------------------------
yields = yEster_cal;
% Example:
% yields = readmatrix('fama_bliss_yields.csv');   % T x N, decimal yields

maturities = mat_Estr_cal(1,:);
dt = 1 / 360;                   % monthly sampling in years
modelType = 'correlated';       % 'correlated' or 'independent'

opts = struct();
opts.maxIter = 20000;
opts.maxFunEvals = 20000;
opts.display = 'iter';
opts.numAdjustmentNodes = 64;
opts.initialParams = [];

if isempty(yields)
    error('Set `yields` to a T x N matrix of zero-coupon yields in decimals before running.');
end

result = estimateAFNS(yields, maturities, dt, modelType, opts);

disp('Estimated AFNS parameters:');
disp(result.params);
disp(['Log likelihood: ', num2str(result.logLik, '%.6f')]);


figure(1);
subplot(2,2,1)
k = 4;
plot(t_y(2:end), yields(2:end,k), '.', 'DisplayName', 'Observed yield');
hold on;
plot(t_y(2:end), result.fittedYields(2:end,k), 'DisplayName', 'AFNS fitted yield');
xlabel('t');
hold off
ylabel('Yield');
legend('Location', 'best');
title('1m yield')
grid on;
subplot(2,2,2)
k = 8;
plot(t_y(2:end), yields(2:end,k), '.', 'DisplayName', 'Observed mean');
hold on;
plot(t_y(2:end), result.fittedYields(2:end,k), 'DisplayName', 'AFNS fitted mean');
xlabel('t');
hold off
ylabel('Yield');
grid on;
title('1y yield')
subplot(2,2,3)
k = 10;
plot(t_y(2:end), yields(2:end,k), '.', 'DisplayName', 'Observed mean');
hold on;
plot(t_y(2:end), result.fittedYields(2:end,k), 'DisplayName', 'AFNS fitted mean');
xlabel('t');
hold off
ylabel('Yield');
grid on;
title('3y yield')
subplot(2,2,4)
k = 17;
plot(t_y(2:end), yields(2:end,k), '.', 'DisplayName', 'Observed mean');
hold on;
plot(t_y(2:end), result.fittedYields(2:end,k), 'DisplayName', 'AFNS fitted mean');
xlabel('t');
hold off
ylabel('Yield');
grid on;
title('10y yield')

function result = estimateAFNS(yields, maturities, dt, modelType, opts)
    if nargin < 5
        opts = struct();
    end
    if ~isfield(opts, 'maxIter'), opts.maxIter = 4000; end
    if ~isfield(opts, 'maxFunEvals'), opts.maxFunEvals = 20000; end
    if ~isfield(opts, 'display'), opts.display = 'iter'; end
    if ~isfield(opts, 'numAdjustmentNodes'), opts.numAdjustmentNodes = 64; end
    if ~isfield(opts, 'initialParams'), opts.initialParams = []; end
    if ~ismember(modelType, {'correlated', 'independent'})
        error('modelType must be ''correlated'' or ''independent''.');
    end

    [numObs, numMat] = size(yields);
    tau = maturities(:);
    B = measurementLoadings(tau, 0.7);

    if ~isempty(opts.initialParams)
        p0 = opts.initialParams(:);
    else
        p0 = initialGuess(yields, modelType);
    end

    objective = @(p) negLogLikAFNS(p, yields, tau, dt, modelType, opts.numAdjustmentNodes);

    optimOpts = optimset( ...
        'Display', opts.display, ...
        'MaxIter', opts.maxIter, ...
        'MaxFunEvals', opts.maxFunEvals, ...
        'TolX', 1e-6, ...
        'TolFun', 1e-6);

    [pHat, fval, exitflag, output] = fminsearch(objective, p0, optimOpts);
    [~, details] = objective(pHat);

    result = struct();
    result.exitflag = exitflag;
    result.output = output;
    result.logLik = -fval;
    result.params = details.model;
    result.filteredStates = details.filteredStates;
    result.predictedStates = details.predictedStates;
    result.filteredCovariances = details.filteredCovariances;
    result.innovations = details.innovations;
    result.fittedYields = details.fittedYields;
    result.residuals = yields - details.fittedYields;
    result.maturities = maturities(:);
    result.dt = dt;
    result.modelType = modelType;
    result.initialMeasurementLoadings = B;
end

function [nll, details] = negLogLikAFNS(p, yields, tau, dt, modelType, numAdjustmentNodes)
    details = struct();

    try
        model = unpackParams(p, numel(tau), modelType);
    catch
        nll = penaltyValue();
        return;
    end

    if ~(model.lambda > 0)
        nll = penaltyValue();
        return;
    end

    eigK = eig(model.KP);
    if any(real(eigK) <= 1e-8)
        nll = penaltyValue();
        return;
    end

    [Phi, Q] = discretizeAFNS(model.KP, model.Sigma, dt);
    if any(~isfinite(Phi(:))) || any(~isfinite(Q(:)))
        nll = penaltyValue();
        return;
    end

    B = measurementLoadings(tau, model.lambda);
    A = afnsYieldAdjustment(tau, model.lambda, model.Sigma, numAdjustmentNodes);
    %H = diag(model.measStd.^ 2);
    H = model.measStd^2*eye(size(yields,2));
    c = (eye(3) - Phi) * model.thetaP;

    [P0ok, P0] = stationaryCovariance(model.KP, model.Sigma);
    if ~P0ok
        nll = penaltyValue();
        return;
    end

    x = model.thetaP;
    P = P0;
    [numObs, numMat] = size(yields);
    logLik = 0;

    filteredStates = zeros(numObs, 3);
    predictedStates = zeros(numObs, 3);
    fittedYields = zeros(numObs, numMat);
    innovations = zeros(numObs, numMat);
    filteredCovariances = zeros(3, 3, numObs);

    for t = 1:numObs
        xPred = c + Phi * x;
        PPred = Phi * P * Phi' + Q;
        PPred = 0.5 * (PPred + PPred');

        yHat = (A + B * xPred).';
        v = yields(t, :) - yHat;
        F = B * PPred * B' + H;
        F = 0.5 * (F + F');

        [cholF, cholFlag] = chol(F, 'lower');
        if cholFlag ~= 0
            nll = penaltyValue();
            return;
        end

        logDetF = 2 * sum(log(diag(cholF)));
        z = cholF \ v.';
        quadForm = z' * z;
        K = (PPred * B') / F;

        x = xPred + K * v.';
        P = PPred - K * B * PPred;
        P = 0.5 * (P + P');

        logLik = logLik - 0.5 * (numMat * log(2 * pi) + logDetF + quadForm);

        predictedStates(t, :) = xPred.';
        filteredStates(t, :) = x.';
        fittedYields(t, :) = yHat;
        innovations(t, :) = v;
        filteredCovariances(:, :, t) = P;
    end

    nll = -logLik;
    details.model = model;
    details.filteredStates = filteredStates;
    details.predictedStates = predictedStates;
    details.filteredCovariances = filteredCovariances;
    details.innovations = innovations;
    details.fittedYields = fittedYields;
end

function model = unpackParams(p, numMat, modelType)
    p = p(:);
    idx = 1;

    model.lambda = exp(p(idx));
    idx = idx + 1;

    if strcmp(modelType, 'independent')
        kpDiag = exp(p(idx:idx + 2));
        model.KP = diag(kpDiag);
        idx = idx + 3;
    else
        kpDiag = exp(p(idx:idx + 2));
        idx = idx + 3;
        kpOff = p(idx:idx + 5);
        idx = idx + 6;
        model.KP = [kpDiag(1), kpOff(1), kpOff(2); ...
                    kpOff(3), kpDiag(2), kpOff(4); ...
                    kpOff(5), kpOff(6), kpDiag(3)];
    end

    model.thetaP = p(idx:idx + 2);
    idx = idx + 3;

    if strcmp(modelType, 'independent')
        sigmaDiag = exp(p(idx:idx + 2));
        idx = idx + 3;
        model.Sigma = diag(sigmaDiag);
    else
        s11 = exp(p(idx)); idx = idx + 1;
        s21 = p(idx);      idx = idx + 1;
        s22 = exp(p(idx)); idx = idx + 1;
        s31 = p(idx);      idx = idx + 1;
        s32 = p(idx);      idx = idx + 1;
        s33 = exp(p(idx)); idx = idx + 1;
        model.Sigma = [s11, 0,   0; ...
                       s21, s22, 0; ...
                       s31, s32, s33];
    end
    model.measStd = exp(p(idx));
%    model.measStd = exp(p(idx:idx + numMat - 1));
end

function p0 = initialGuess(yields, modelType)
    yBar = mean(yields, 1);
    level0 = yBar(end);
    slope0 = yBar(1) - yBar(end);
    curve0 = 2 * yBar(round(end / 2)) - yBar(1) - yBar(end);
    theta0 = [level0; slope0; curve0];

    lambda0 = log(0.7);
    %meas0 = log(1e-3 * ones(size(yields, 2), 1));
    meas0 = log(1e-3);

    if strcmp(modelType, 'independent')
        kp0 = log([0.08; 0.20; 1.20]);
        sigma0 = log([0.01; 0.01; 0.01]);
        p0 = [lambda0; kp0; theta0; sigma0; meas0];
    else
        kp0 = [log(0.20); log(0.40); log(1.20); 0; 0; 0; 0; 0; 0];
        sigma0 = [log(0.01); 0; log(0.01); 0; 0; log(0.01)];
        p0 = [lambda0; kp0; theta0; sigma0; meas0];
    end
end

function B = measurementLoadings(tau, lambda)
    x = lambda * tau;
    slope = safeNS1(x);
    curvature = slope - exp(-x);
    B = [ones(numel(tau), 1), slope, curvature];
end

function adj = afnsYieldAdjustment(tau, lambda, Sigma, numNodes)
    SigmaSigmaT = Sigma * Sigma.';
    adj = zeros(numel(tau), 1);

    for i = 1:numel(tau)
        Ti = tau(i);
        nodes = linspace(0, Ti, numNodes);
        B1 = -(Ti - nodes);
        x = lambda * (Ti - nodes);
        B2 = -safeNS1(x) .* (Ti - nodes);
        B3 = (Ti - nodes) .* exp(-x) - safeNS1(x) .* (Ti - nodes);

        integrand = zeros(size(nodes));
        for k = 1:numel(nodes)
            b = [B1(k), B2(k), B3(k)];
            integrand(k) = b * SigmaSigmaT * b.';
        end

        bigA = 0.5 * trapz(nodes, integrand);
        adj(i) = -bigA / Ti;
    end
end

function y = safeNS1(x)
    y = zeros(size(x));
    small = abs(x) < 1e-8;
    y(small) = 1 - x(small) / 2 + x(small).^2 / 6;
    y(~small) = (1 - exp(-x(~small))) ./ x(~small);
end

function [Phi, Q] = discretizeAFNS(KP, Sigma, dt)
    Phi = expm(-KP * dt);
    A = [-KP, Sigma * Sigma.'; zeros(3), KP.'] * dt;
    M = expm(A);
    Q = M(1:3, 4:6) * Phi.';
    Q = 0.5 * (Q + Q.');
end

function [ok, P] = stationaryCovariance(KP, Sigma)
    M = kron(eye(3), KP) + kron(KP, eye(3));
    rhs = reshape(Sigma * Sigma.', [], 1);

    if rcond(M) < 1e-12
        ok = false;
        P = nan(3);
        return;
    end

    P = reshape(M \ rhs, 3, 3);
    P = 0.5 * (P + P.');
    eigP = eig(P);
    ok = all(isfinite(P(:))) && all(real(eigP) > -1e-10);
end

function val = penaltyValue()
    val = 1e12;
end
