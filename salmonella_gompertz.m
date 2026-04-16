%% Salmonella Gompertz Growth Model - Sinusoidal Temperature
% Primary model (differential): d(logN)/dt = mu*K*C*exp(-K)
% where K = exp(-mu*(t-M)), computed directly at each time step
% Secondary model: mu = a*(T-Tmin)^2*(1-exp(b*(T-Tmax)))
%%
%% Preprocessing
clear
close all; 
format compact
global tTemp
figDir = fullfile(fileparts(mfilename('fullpath')), 'figures');
if ~exist(figDir, 'dir'), mkdir(figDir); end
%%
%% Read in data
% --- Growth data (logN) ---
growthData = readmatrix('Salmonella sin growth.xlsx');
% readmatrix skips text headers automatically
% Col1=time(hr), Col2=CFU/mL rep1, Col3=CFU/mL rep2
% times 0,0.5,1,2 hr have only 1 replicate; others have 2

rawTime = growthData(:,1);
rawCFU1 = growthData(:,2);
rawCFU2 = growthData(:,3);

% Build combined vectors: duplicate times for replicates
x = [];      % time vector for logN data
yobs = [];   % logN observed
for i = 1:length(rawTime)
    if ~isnan(rawCFU1(i)) && rawCFU1(i) > 0
        x    = [x; rawTime(i)];
        yobs = [yobs; log10(rawCFU1(i))];
    end
    if ~isnan(rawCFU2(i)) && rawCFU2(i) > 0
        x    = [x; rawTime(i)];
        yobs = [yobs; log10(rawCFU2(i))];
    end
end
n = length(x);
fprintf('n = %d data points\n', n);

% --- Temperature data ---
tempData = readmatrix('Salmonella sin growth Temps.xlsx');
tTemp_min = tempData(:,1);       % time in minutes
tTemp_hr  = tTemp_min / 60;      % convert to hours
Tobs      = tempData(:,2);       % temperature degC
tTemp     = [tTemp_hr, Tobs];    % global: [time(hr), T(degC)]
fprintf('Temp time range: [%.2f, %.2f] hr\n', min(tTemp_hr), max(tTemp_hr));
%%
%% Initial parameter guesses -- ALL 7 parameters
% beta_all = [A, C, M, a, b, Tmin, Tmax]
A    = log10(400);  % log(400 cfu/ml)
C    = 11;          % log cfu/ml
M    = 7.5;         % hr
a    = 0.000338;    % degC^-2 hr^-1
b    = 0.275;       % degC^-1
Tmin = 6;           % degC
Tmax = 46.3;        % degC

beta_all = [A, C, M, a, b, Tmin, Tmax];
pnames_all = {'A','C','M','a','b','T_{min}','T_{max}'};
p_all = length(beta_all);

paramColorMap = containers.Map( ...
    {'A','C','M','a','b','T_{min}','T_{max}'}, ...
    {[1 0 0], [0 0.6 0], [0 0 1], [0 0.7 0.7], [0.85 0.65 0], [0.7 0 0.7], [0.3 0.3 0.3]});
Ypred_color = [0 0 0];
%%
%% define the function for forward problem (all 7 params)
fnameFOR_all = @gompertzFOR_7;
%%
%% plot ypred with initial parameter guesses on data
% to make sure guesses are reasonable
xs = linspace(min(x), max(x), 200)';
ns = length(xs);
ypredInit = fnameFOR_all(beta_all, xs);
figure
set(gca, 'fontsize',16,'fontweight','bold');
plot(xs, ypredInit, '-k', 'LineWidth', 2);
hold on
plot(x, yobs, 'or', 'MarkerSize', 8, 'MarkerFaceColor', 'r');
hold off
xlabel('time (hr)'); ylabel('log_{10}N (log cfu/mL)');
title('Initial guess vs observed data');
legend('Predicted (initial guess)','Observed','location','best');
saveas(gcf, fullfile(figDir, 'fig01_initial_guess_vs_data.png'));
%%
%% X' = scaled sensitivity coefficients (ALL 7 parameters)
% This is a forward problem with known approximate parameters
Xp_all = SSC_V3(beta_all, xs, fnameFOR_all);
%%
%% plot X' for all 7 parameters
figure
hold on
set(gca, 'fontsize',18,'fontweight','bold');
for i = 1:p_all
    h2(i) = plot(xs(1:ns), Xp_all(1:ns,i), '-', 'color', paramColorMap(pnames_all{i}), 'LineWidth', 3);
end
ypred = fnameFOR_all(beta_all, xs);
h2(p_all+1) = plot(xs, ypred, '--', 'color', Ypred_color, 'LineWidth', 4);
legStr = cell(1, p_all+1);
for i = 1:p_all
    legStr{i} = [pnames_all{i}, '*\partialY/\partial(', pnames_all{i}, ')'];
end
legStr{p_all+1} = 'Y';
legend(h2, legStr, 'location', 'best');
xlabel('time (hr)');
ylabel('SSC  \beta_i \partial Y/\partial\beta_i  (Y units)');
title('SSC using initial guesses -- ALL 7 parameters');
grid on
saveas(gcf, fullfile(figDir, 'fig02_SSC_7params.png'));

%% SSC ratio analysis -- skipped (collinearity now checked via correlation matrix)

% Print max |SSC| for each parameter to help decide which to fix
fprintf('\n--- SSC Analysis (all 7 params) ---\n');
maxSSC_all = max(abs(Xp_all));
for i = 1:p_all
    fprintf('  Max |SSC(%s)| = %.4f\n', pnames_all{i}, maxSSC_all(i));
end
[~, sscOrder] = sort(maxSSC_all, 'descend');
fprintf('\n4b. Parameters ranked by estimation accuracy (highest SSC = most accurately estimated):\n');
for k = 1:p_all
    fprintf('  %d. %s  (Max|SSC|=%.4f)\n', k, pnames_all{sscOrder(k)}, maxSSC_all(sscOrder(k)));
end
fprintf('\nBased on SSC analysis, Tmin (0.83) and b (0.88) have lowest SSC => fix them.\n');
fprintf('Fix Tmin = %.1f degC, b = %.6f at literature values.\n', Tmin, b);
fprintf('Remaining parameters to estimate: A, C, M, a, Tmax (p=5)\n');
%%
%% nlinfit options (shared across all rounds and bootstrap)
opts = statset('nlinfit');
opts.MaxIter = 1000;
%%
%% ======== Round 0: Estimate ALL 7 parameters (for comparison) ========
fnameINV_all = @gompertzFOR_7;  % same function works for inverse
fprintf('\n===== Round 0: Attempting 7-parameter estimation =====\n');
try
    [beta_r0, resids_r0, J_r0, COVB_r0, mse_r0] = nlinfit(x, yobs, fnameINV_all, beta_all, opts);
    fprintf('beta_r0 =\n'); disp(beta_r0);
    rmse_r0 = sqrt(mse_r0);
    SSres_r0 = resids_r0'*resids_r0;
    SStot_r0 = sum((yobs - mean(yobs)).^2);
    Rsq_r0 = 1 - SSres_r0/SStot_r0;
    fprintf('R^2     = %.4f\n', Rsq_r0);
    fprintf('RMSE    = %.4f\n', rmse_r0);
    fprintf('cond(J) = %.4e\n', cond(J_r0));
    ci_r0 = nlparci(beta_r0, resids_r0, J_r0);
    fprintf('95%% CI:\n'); 
    for i = 1:p_all
        fprintf('  %s: [%.6g, %.6g]\n', pnames_all{i}, ci_r0(i,1), ci_r0(i,2));
    end
    
    % Plot 7-param fit
    figure
    hold on
    set(gca, 'fontsize',14,'fontweight','bold');
    ypredp_r0 = fnameFOR_all(beta_r0, xs);
    plot(xs, ypredp_r0, '-b', 'LineWidth', 2);
    plot(x, yobs, 'sr', 'MarkerSize', 8, 'MarkerFaceColor', 'r');
    xlabel('time (hr)'); ylabel('log_{10}N (log cfu/mL)');
    title(sprintf('7-param fit: R^2=%.4f, RMSE=%.4f, cond(J)=%.2e', Rsq_r0, rmse_r0, cond(J_r0)));
    legend('7-param predicted', 'Observed', 'location', 'best');
    grid on
    saveas(gcf, fullfile(figDir, 'fig00_7param_fit.png'));
catch ME
    fprintf('7-parameter estimation failed: %s\n', ME.message);
end
fprintf('Note: 7-param results shown for comparison only.\n');
%%
%% ======== Round 1: Fix Tmin, b -- estimate 5 parameters [A,C,M,a,Tmax] ========
Tmin_fixed = Tmin;
b_fixed = b;

beta0_r1(1) = A;       % A
beta0_r1(2) = C;       % C
beta0_r1(3) = M;       % M
beta0_r1(4) = a;       % a
beta0_r1(5) = Tmax;    % Tmax
p_r1 = length(beta0_r1);
pnames_r1 = {'A','C','M','a','T_{max}'};
%%
%% define functions for Round 1 (5 params, Tmin & b fixed)
fnameFOR_r1 = @(beta,t) gompertzFOR_5(beta, t, Tmin_fixed, b_fixed);
fnameINV_r1 = @(beta,t) gompertzINV_5(beta, t, Tmin_fixed, b_fixed);
%%
%% plot ypred with 5-param initial guesses on data
ypredInit = fnameFOR_r1(beta0_r1, xs);
figure
plot(xs, ypredInit, '-k', x, yobs, 'or');
xlabel('time (hr)'); ylabel('log_{10}N (log cfu/mL)');
title('Initial guess (5 params, Tmin & b fixed) vs observed data');
legend('Predicted (initial guess)','Observed','location','best');
% saveas(gcf, fullfile(figDir, 'fig03_initial_guess_5params.png'));
%%
%% X' = SSC for 5 estimable parameters (Round 1)
Xp_r1 = SSC_V3(beta0_r1, xs, fnameFOR_r1);
%%
%% plot X' for 5 parameters (Round 1)
figure
hold on
set(gca, 'fontsize',18,'fontweight','bold');
clear h2
for i = 1:p_r1
    h2(i) = plot(xs(1:ns), Xp_r1(1:ns,i), '-', 'color', paramColorMap(pnames_r1{i}), 'LineWidth', 3);
end
ypred = fnameFOR_r1(beta0_r1, xs);
h2(p_r1+1) = plot(xs, ypred, '--', 'color', Ypred_color, 'LineWidth', 4);
legStr_r1 = cell(1, p_r1+1);
for i = 1:p_r1
    legStr_r1{i} = [pnames_r1{i}, '*\partialY/\partial(', pnames_r1{i}, ')'];
end
legStr_r1{p_r1+1} = 'Y';
legend(h2, legStr_r1, 'location', 'best');
xlabel('time (hr)');
ylabel('SSC  \beta_i \partial Y/\partial\beta_i  (Y units)');
title('SSC using initial guesses -- 5 params (Tmin & b fixed)');
grid on
saveas(gcf, fullfile(figDir, 'fig06a_SSC_5params_initial.png'));
%%
%% Round 1 nlinfit -- OLS inverse problem (5 params)
[beta_r1, resids_r1, J_r1, COVB_r1, mse_r1] = nlinfit(x, yobs, fnameINV_r1, beta0_r1, opts);
fprintf('\n===== Round 1 Results (5 params) =====\n');
beta_r1
rmse_r1 = sqrt(mse_r1)
condX_r1 = cond(J_r1)
detXTX_r1 = det(J_r1'*J_r1)
%%
%% Round 1 confidence intervals for parameters
ci_r1 = nlparci(beta_r1, resids_r1, J_r1)
[R_r1, sigma_r1] = corrcov(COVB_r1);
R_r1
sigma_r1
relerr_r1 = sigma_r1./beta_r1'
%%
%% Round 1 SSC using estimated parameters -- check Tmax identifiability
clear Xp_r1_final
Xp_r1_final = SSC_V3(beta_r1, xs, fnameFOR_r1);

fprintf('\n--- Round 1 SSC Analysis (5 params, estimated) ---\n');
for i = 1:p_r1
    fprintf('  Max |SSC(%s)| = %.6f\n', pnames_r1{i}, max(abs(Xp_r1_final(:,i))));
end

figure
clear h2
hold on
set(gca, 'fontsize',14,'fontweight','bold');
for i = 1:p_r1
    h2(i) = plot(xs(1:ns), Xp_r1_final(1:ns,i), '-', 'color', paramColorMap(pnames_r1{i}), 'LineWidth', 3);
end
ypred = fnameFOR_r1(beta_r1, xs);
h2(p_r1+1) = plot(xs, ypred, '--', 'color', Ypred_color, 'LineWidth', 4);
legend(h2, legStr_r1, 'location', 'best');
xlabel('time (hr)');
ylabel('SSC  \beta_i \partial Y/\partial\beta_i  (Y units)');
title('SSC using estimated parameters (5 params)');
grid on
% saveas(gcf, fullfile(figDir, 'fig05_SSC_5params_estimated.png'));

fprintf('\n--- Collinearity Check (5 params) ---\n');
for ii = 1:p_r1
    for jj = ii+1:p_r1
        fprintf('  R(%s, %s) = %.4f', pnames_r1{ii}, pnames_r1{jj}, R_r1(ii,jj));
        if abs(R_r1(ii,jj)) > 0.95
            fprintf('  <-- highly correlated');
        end
        fprintf('\n');
    end
end
fprintf('  cond(J) = %.4e  (should be < 1e6)\n', condX_r1);
fprintf('R(C, Tmax)=%.3f: correlated but still estimable (95%% CI excludes zero).\n', R_r1(2,5));
fprintf('Final estimable parameters: A, C, M, a, Tmax (p=5)\n');
%%
%% ======== Final: Round 1 is the final estimation (5 params) ========
% Reviewer: R(C,Tmax)~0.97 is just enough to allow estimation;
% 95% CI excludes zero even though error is larger.
beta = beta_r1;
resids = resids_r1;
J = J_r1;
COVB = COVB_r1;
mse = mse_r1;
p = p_r1;
pnames = pnames_r1;
fnameFOR = fnameFOR_r1;
fnameINV = fnameINV_r1;
n = length(x);
rmse = sqrt(mse);
condX = cond(J);
detXTX = det(J'*J);

SSres = resids'*resids;
SStot = sum((yobs - mean(yobs)).^2);
Rsq = 1 - SSres/SStot;
Rsq_adj = 1 - (SSres/(n-p)) / (SStot/(n-1));

fprintf('\n===== Final OLS Results (5 params: [A,C,M,a,Tmax]) =====\n');
fprintf('beta =\n'); disp(beta);
fprintf('R^2     = %.4f\n', Rsq);
fprintf('R^2_adj = %.4f\n', Rsq_adj);
fprintf('RMSE    = %.4f\n', rmse);
fprintf('relRMSE = %.4f\n', rmse/range(fnameINV(beta,x)));
fprintf('cond(J) = %.4e\n', condX);
fprintf('det(J''J) = %.4e\n', detXTX);

ci = nlparci(beta, resids, J);
fprintf('ci =\n'); disp(ci);
[R, sigma] = corrcov(COVB);
fprintf('R =\n'); disp(R);
fprintf('sigma =\n'); disp(sigma);
relerr = sigma./beta';
fprintf('relerr =\n'); disp(relerr);
%%
%% Confidence and prediction intervals for the dependent variable
[ypred, delta]   = nlpredci(fnameINV, x, beta, resids, J, 0.05, 'on', 'curve');
[ypred, deltaob] = nlpredci(fnameINV, x, beta, resids, J, 0.05, 'on', 'observation');

yspan = range(ypred)
relrmse = rmse/yspan

% simultaneous confidence bands for regression line
CBu = ypred + delta;
CBl = ypred - delta;
% simultaneous prediction bands for individual points
PBu = ypred + deltaob;
PBl = ypred - deltaob;
%%
%% Output -- ypred and yobs vs. t
figure
hold on
ypredp = fnameFOR(beta, xs); % smooth curve for plotting
set(gca, 'fontsize',20,'fontweight','bold');
h1(1) = plot(xs, ypredp, '-', 'linewidth', 3);
h1(2) = plot(x, yobs, 'square', 'Markerfacecolor', 'r');
xlabel('time (hr)','fontsize',16,'fontweight','bold')
ylabel('log_{10}N (log cfu/mL)','fontsize',16,'fontweight','bold')
%%
%% Output -- CIs and PIs
[xSort, idx] = sort(x);
CBuS = CBu(idx); CBlS = CBl(idx);
PBuS = PBu(idx); PBlS = PBl(idx);

h1(3) = plot(xSort, CBuS, '--g', 'LineWidth', 2);
plot(xSort, CBlS, '--g', 'LineWidth', 2);
h1(4) = plot(xSort, PBuS, '-.c', 'LineWidth', 2);
plot(xSort, PBlS, '-.c', 'LineWidth', 2);
legend(h1, 'ypred', 'yobs', 'CB', 'PB')
title('OLS fit with asymptotic CB and PB')
saveas(gcf, fullfile(figDir, 'fig06_OLS_fit_CB_PB.png'));
%%
%% 7d. Dual Y-axis plot: logN observed + predicted (left) and Temperature (right)
figure
set(gca, 'fontsize',16,'fontweight','bold');
yyaxis left
hd(1) = plot(x, yobs, 'sb', 'MarkerSize', 8, 'MarkerFaceColor', 'b');
hold on
hd(2) = plot(xs, ypredp, '-b', 'LineWidth', 2);
ylabel('log_{10}N (log cfu/mL)');
yyaxis right
hd(3) = plot(tTemp_hr, Tobs, '-r', 'LineWidth', 2);
ylabel('Temperature (°C)');
xlabel('time (hr)');
title('Salmonella growth under sinusoidal temperature');
legend(hd, 'logN observed', 'logN predicted', 'Temperature', 'location', 'best');
grid on
saveas(gcf, fullfile(figDir, 'fig07_dual_axis_logN_temp.png'));
figure
hold on
plot(x, resids, 'square', 'Markerfacecolor', 'b', 'markersize', 10);
YLine = [0 0];
XLine = [0 max(x)];
plot(XLine, YLine, 'R');
ylabel('Observed y - Predicted y','fontsize',16,'fontweight','bold')
xlabel('Time (hr)','fontsize',16,'fontweight','bold')
grid on
set(gca, 'fontsize',14,'fontweight','bold');
meanres = mean(resids)
hold off
saveas(gcf, fullfile(figDir, 'fig08_residual_scatter.png'));
%%
%% number of runs (handles replicates)
x3 = x;
xResids = [x3 resids];
xResidsSort = sortrows(xResids);
xS = xResidsSort(:,1);
residsSort = xResidsSort(:,2);

count = 0;
countRep = 1;
countNeg = 0; countPos = 0;
resSign(countRep) = sign(residsSort(1));
for i = 2:n
    if xS(i)==xS(i-1)
        countRep = countRep+1;
        resSign(countRep) = sign(residsSort(i));
    else
        for j = 1:countRep
            if resSign(j) < 0
                countNeg = countNeg+1;
            elseif resSign(j) > 0
                countPos = countPos+1;
            end
        end
        count = count+min(countNeg, countPos);
        rescross = residsSort(i)*residsSort(i-1);
        resSign(j+1) = sign(rescross);
        if resSign(j+1) < 0
            count = count+1;
        end
        clear resSign
        countRep = 1;
        resSign(countRep) = sign(residsSort(i));
        countNeg = 0; countPos = 0;
    end
end
fprintf('number of runs = %5.2f\n', count);
minrun = (n+1)/2;
fprintf('Minimum required number of runs = %5.2f\n', minrun);
%%
%% residuals histogram
figure
h = histogram(resids);
set(gca, 'fontsize',14,'fontweight','bold');
title('Residual Histogram');
xlabel('Y_{observed} - Y_{predicted}','fontsize',16,'fontweight','bold')
ylabel('Frequency','fontsize',16,'fontweight','bold')
saveas(gcf, fullfile(figDir, 'fig09_residual_histogram.png'));
clear Xp ypred
Xp = SSC_V3(beta, xs, fnameFOR);
%%
%% plot X' final (5 params, after nlinfit)
figure
hold on
set(gca, 'fontsize',14,'fontweight','bold');
clear h2
legStr_final = cell(1, p+1);
for i = 1:p
    h2(i) = plot(xs(1:ns), Xp(1:ns,i), '-', 'color', paramColorMap(pnames{i}), 'LineWidth', 3);
    legStr_final{i} = [pnames{i}, '*\partialY/\partial(', pnames{i}, ')'];
end
ypred = fnameFOR(beta, xs);
h2(p+1) = plot(xs, ypred, '--', 'color', Ypred_color, 'LineWidth', 4);
legStr_final{p+1} = 'Y';
legend(h2, legStr_final, 'location', 'best');
xlabel('time (hr)');
ylabel('SSC  \beta_i \partial Y/\partial\beta_i  (Y units)');
maxSSC_final = max(abs(Xp));
title('Final SSC — after nlinfit (5 params)');
grid on
saveas(gcf, fullfile(figDir, 'fig10_SSC_5params_final.png'));
%%
%% Optimal Experimental Design -- Beck & Arnold (1977) Eq. 8.3.5
% C_ij(t) = (1/t) * cumtrapz(t, X'_i * X'_j)
% delta(t) = det(C(t)),  Cii = diag(C(t))
% Follows optexp_Fig8_10v2.m exactly
tOED = linspace(0, max(tTemp_hr), 300)';
mOED = length(tOED);

Xp_oed = SSC_V3(beta, tOED, fnameFOR);

Coed = cell(p, p);
for i = 1:p
    for j = 1:p
        intgrnd = Xp_oed(:,i) .* Xp_oed(:,j);
        Coed{i,j} = (1./tOED) .* cumtrapz(tOED, intgrnd);
    end
end

CC = zeros(p, p, mOED);
for i = 1:p
    for j = 1:p
        CC(i,j,:) = Coed{i,j};
    end
end

CC(:,:,1) = 0;
CC(1,1,1) = 1;  % A (param 1) is the initial value => C11(0) = 1

delta_oed = zeros(1, mOED);
delta_oed(1) = 0;
for k = 2:mOED
    delta_oed(k) = det(CC(:,:,k));
end

Cp = zeros(p, mOED);
for i = 1:p
    Cp(i,:) = CC(i,i,:);
end

figure
subplot(2,1,1)
plot(tOED, delta_oed, '-b', 'LineWidth', 2);
set(gca, 'fontsize',14,'fontweight','bold');
xlabel('time (hr)');
ylabel('\Delta = det(C)');
title('\Delta Criterion  (Beck & Arnold Eq. 8.3.5)');
grid on

subplot(2,1,2)
hold on
set(gca, 'fontsize',14,'fontweight','bold');
for i = 1:p
    plot(tOED, Cp(i,:), '-', 'color', paramColorMap(pnames{i}), 'LineWidth', 2);
end
xlabel('time (hr)');
ylabel('C_{ii}');
title('Diagonal elements of C');
legend(pnames, 'Location', 'best');
grid on
hold off
saveas(gcf, fullfile(figDir, 'fig11_OED_delta_Cii.png'));
nBoot = 600;
betaBoot = zeros(nBoot, p);
ypredBoot = zeros(nBoot, ns);

ypredObs = fnameINV(beta, x); % predicted at observation times

% Bootstrap uses relaxed tolerances: statistical resampling does not need
% machine-precision convergence. Looser TolFun/TolX resolves most
% IterationLimitExceeded warnings that arise from the high C--Tmax and
% A--M parameter correlations distorting the optimization landscape for
% perturbed bootstrap datasets.
optsBoot = statset('nlinfit');
optsBoot.MaxIter = 1000;
optsBoot.TolFun  = 1e-6;   % relaxed from default 1e-8
optsBoot.TolX    = 1e-6;   % relaxed from default 1e-8

rng(42);
nFailed = 0;
warning('off', 'stats:nlinfit:IllConditionedJacobian');
warning('off', 'stats:nlinfit:ModelConstantWRTParam');
warning('off', 'MATLAB:rankDeficientMatrix');
warning('off', 'MATLAB:nearlySingularMatrix');
for ib = 1:nBoot
    bootIdx = randi(n, n, 1);
    bootResids = resids(bootIdx);
    yBoot = ypredObs + bootResids;
    try
        [betaB,~,~,~,~] = nlinfit(x, yBoot, fnameINV, beta, optsBoot);
        betaBoot(ib,:) = betaB;
        ypredBoot(ib,:) = fnameFOR(betaB, xs);
    catch
        betaBoot(ib,:) = NaN;
        ypredBoot(ib,:) = NaN;
        nFailed = nFailed + 1;
    end
    if mod(ib,100)==0
        fprintf('Bootstrap iteration %d/%d\n', ib, nBoot);
    end
end
warning('on', 'stats:nlinfit:IllConditionedJacobian');
warning('on', 'stats:nlinfit:ModelConstantWRTParam');
warning('on', 'MATLAB:rankDeficientMatrix');
warning('on', 'MATLAB:nearlySingularMatrix');
if nFailed > 0
    fprintf('Bootstrap: %d/%d iterations failed (caught as NaN)\n', nFailed, nBoot);
end

validIdx = ~any(isnan(betaBoot),2);
betaBoot = betaBoot(validIdx,:);
ypredBoot = ypredBoot(validIdx,:);

% Remove outliers: discard iterations where any parameter deviates > 5x from estimate
outlierIdx = false(size(betaBoot,1), 1);
for i = 1:p
    ratio = abs(betaBoot(:,i) ./ beta(i));
    outlierIdx = outlierIdx | (ratio > 5) | (ratio < 0.2);
end
if any(outlierIdx)
    fprintf('Bootstrap: removed %d outlier iterations\n', sum(outlierIdx));
    betaBoot = betaBoot(~outlierIdx,:);
    ypredBoot = ypredBoot(~outlierIdx,:);
end
fprintf('Valid bootstrap iterations: %d/%d\n', size(betaBoot,1), nBoot);

fprintf('\n===== Bootstrap 95%% CI =====\n');
for i = 1:p
    bootCI = prctile(betaBoot(:,i), [2.5 97.5]);
    fprintf('  %s: [%.6g, %.6g]\n', pnames{i}, bootCI(1), bootCI(2));
end

% Bootstrap CB and PB
bootCBu = prctile(ypredBoot, 97.5, 1)';
bootCBl = prctile(ypredBoot, 2.5, 1)';
bootPBu = bootCBu + 1.96*rmse;
bootPBl = bootCBl - 1.96*rmse;

figure
hold on
set(gca, 'fontsize',16,'fontweight','bold');
hb(1) = plot(xs, fnameFOR(beta,xs), '-b', 'LineWidth', 2);
hb(2) = plot(x, yobs, 'sr', 'MarkerSize', 8, 'MarkerFaceColor', 'r');
hb(3) = plot(xs, bootCBu, '--g', 'LineWidth', 2);
plot(xs, bootCBl, '--g', 'LineWidth', 2);
hb(4) = plot(xs, bootPBu, '-.c', 'LineWidth', 2);
plot(xs, bootPBl, '-.c', 'LineWidth', 2);
xlabel('time (hr)');
ylabel('log_{10}N (log cfu/mL)');
title('Bootstrap CB and PB');
legend(hb, 'ypred', 'yobs', 'Boot CB', 'Boot PB', 'location', 'best');
grid on
hold off
saveas(gcf, fullfile(figDir, 'fig12_bootstrap_CB_PB.png'));
fprintf('\n===== Asymptotic vs Bootstrap Bands =====\n');
fprintf('Asymptotic CB width (avg): %.4f\n', mean(2*delta));
fprintf('Bootstrap  CB width (avg): %.4f\n', mean(bootCBu - bootCBl));
fprintf('Asymptotic PB width (avg): %.4f\n', mean(2*deltaob));
fprintf('Bootstrap  PB width (avg): %.4f\n', mean(bootPBu - bootPBl));

%%
%% ==================== FUNCTIONS ====================

%% Forward model -- ALL 7 parameters (for initial SSC analysis)
function logN = gompertzFOR_7(beta, t)
    A=beta(1); C=beta(2); M=beta(3);
    a=beta(4); b=beta(5); Tmin=beta(6); Tmax=beta(7);
    tAll = t(:);
    [tUniq, ~, ic] = unique(tAll);
    y0 = A;
    [~,Y] = ode45(@(t,y) gompODE(t,y,C,M,a,b,Tmin,Tmax), tUniq, y0);
    logN = Y(ic, 1);
end

%% Forward model -- 5 parameters (Tmin & b fixed) -- Round 1
function logN = gompertzFOR_5(beta, t, Tmin, b)
    A=beta(1); C=beta(2); M=beta(3);
    a=beta(4); Tmax=beta(5);
    tAll = t(:);
    [tUniq, ~, ic] = unique(tAll);
    y0 = A;
    [~,Y] = ode45(@(t,y) gompODE(t,y,C,M,a,b,Tmin,Tmax), tUniq, y0);
    logN = Y(ic, 1);
end

%% Inverse model -- 5 parameters (Tmin & b fixed) -- Round 1
function logN = gompertzINV_5(beta, t, Tmin, b)
    A=beta(1); C=beta(2); M=beta(3);
    a=beta(4); Tmax=beta(5);
    tAll = t(:);
    [tUniq, ~, ic] = unique(tAll);
    y0 = A;
    [~,Y] = ode45(@(t,y) gompODE(t,y,C,M,a,b,Tmin,Tmax), tUniq, y0);
    logN = Y(ic, 1);
end

%% Forward model -- 4 parameters (Tmin, b, Tmax fixed) -- Round 2
function logN = gompertzFOR_4a(beta, t, Tmin, b, Tmax)
    A=beta(1); C=beta(2); M=beta(3);
    a=beta(4);
    tAll = t(:);
    [tUniq, ~, ic] = unique(tAll);
    y0 = A;
    [~,Y] = ode45(@(t,y) gompODE(t,y,C,M,a,b,Tmin,Tmax), tUniq, y0);
    logN = Y(ic, 1);
end

%% Inverse model -- 4 parameters (Tmin, b, Tmax fixed) -- Round 2
function logN = gompertzINV_4a(beta, t, Tmin, b, Tmax)
    A=beta(1); C=beta(2); M=beta(3);
    a=beta(4);
    tAll = t(:);
    [tUniq, ~, ic] = unique(tAll);
    y0 = A;
    [~,Y] = ode45(@(t,y) gompODE(t,y,C,M,a,b,Tmin,Tmax), tUniq, y0);
    logN = Y(ic, 1);
end

%% ODE: d(logN)/dt = mu*K*C*exp(-K)
% K = exp(-mu*(t-M)) computed directly, NOT a separate state variable
function dy = gompODE(t, y, C, M, a, b, Tmin, Tmax)
    global tTemp
    T = interp1(tTemp(:,1), tTemp(:,2), t, 'linear', 'extrap');
    mu = secModel(T, a, b, Tmin, Tmax);
    K = exp(-mu*(t - M));
    dy = mu*K*C*exp(-K);
end

%% Secondary model: mu = a*(T-Tmin)^2*(1-exp(b*(T-Tmax)))
function mu = secModel(T, a, b, Tmin, Tmax)

    mu = a*(T-Tmin).^2.*(1-exp(b*(T-Tmax)));

end

%% SSC_V3 -- scaled sensitivity coefficients via forward difference
function Xp = SSC_V3(beta, x, yfunc)
    d = 0.001;
    ypred = yfunc(beta, x);
    for i = 1:length(beta)
        betain = beta;
        betain(i) = beta(i)*(1+d);
        yhat = yfunc(betain, x);
        Xp(:,i) = (yhat - ypred)/d;
    end
end

%% JACOB_CD -- unscaled Jacobian via central difference (matches nlinfit approach)
function J = JACOB_CD(beta, x, yfunc)
    d = 0.001;
    p = length(beta);
    n = length(x);
    J = zeros(n, p);
    for i = 1:p
        betaP = beta; betaP(i) = beta(i)*(1+d);
        betaM = beta; betaM(i) = beta(i)*(1-d);
        yP = yfunc(betaP, x);
        yM = yfunc(betaM, x);
        J(:,i) = (yP - yM) / (2*d*beta(i));
    end
end
