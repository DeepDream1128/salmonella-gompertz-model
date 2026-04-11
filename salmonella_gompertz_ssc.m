%[text] # Salmonella Gompertz -- Scaled Sensitivity Coefficient (SSC) Analysis
% This script focuses on SSC results for the Gompertz growth model
% under sinusoidal temperature conditions.
%%
%[text] ## Preprocessing
s = settings;
s.matlab.fonts.editor.code.Size.TemporaryValue = '14pt';
clear
close all;
format compact
global tTemp
%%
%[text] ## Read in data
% --- Growth data (logN) ---
growthData = readmatrix('Salmonella sin growth.xlsx');
rawTime = growthData(:,1);
rawCFU1 = growthData(:,2);
rawCFU2 = growthData(:,3);

x = [];
yobs = [];
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

% --- Temperature data ---
tempData = readmatrix('Salmonella sin growth Temps.xlsx');
tTemp_min = tempData(:,1);
tTemp_hr  = tTemp_min / 60;
Tobs      = tempData(:,2);
tTemp     = [tTemp_hr, Tobs];
%%
%[text] ## Initial parameter guesses -- ALL 7 parameters
A    = log10(400);
C    = 11;
M    = 7.5;
a    = 0.000338;
b    = 0.275;
Tmin = 6;
Tmax = 46.3;

beta_all = [A, C, M, a, b, Tmin, Tmax];
pnames_all = {'A','C','M','a','b','T_{min}','T_{max}'};
p_all = length(beta_all);
%%
%[text] ## define the function that will be used for the forward problem (all 7 params)
fnameFOR_all = @gompertzFOR_7;
%%
%[text] ## X' = SSC using initial guesses -- ALL 7 parameters
xs = linspace(min(x), max(x), 200)';
ns = length(xs);
Xp_all = SSC_V3(beta_all, xs, fnameFOR_all);
%%
%[text] ## plot X' for all 7 parameters
cmap = ['r' 'g' 'b' 'c' 'y' 'm' 'k']';
figure
hold on
set(gca, 'fontsize',18,'fontweight','bold');
for i = 1:p_all
    h2(i) = plot(xs(1:ns), Xp_all(1:ns,i), '-', 'color', cmap(i,:), 'LineWidth', 3);
end
ypred = fnameFOR_all(beta_all, xs);
h2(p_all+1) = plot(xs, ypred, '--k', 'LineWidth', 4);
legStr_all = cell(1, p_all+1);
for i = 1:p_all
    legStr_all{i} = [pnames_all{i}, '*\partialY/\partial(', pnames_all{i}, ')'];
end
legStr_all{p_all+1} = 'Y';
legend(h2, legStr_all, 'location', 'best');
xlabel('time (hr)');
ylabel('scaled sensitivity coefficient \beta_i*\partial(Y)/\partial(\beta_i), Y units');
title('SSC using initial guesses -- ALL 7 parameters');
grid on
%%
%[text] ## SSC magnitude analysis -- decide which parameters to fix
fprintf('\n========== SSC Analysis (all 7 params) ==========\n');
for i = 1:p_all
    fprintf('  Max |SSC(%s)| = %.6f\n', pnames_all{i}, max(abs(Xp_all(:,i))));
end
fprintf('\nTmin and Tmax SSCs are near zero => cannot be estimated.\n');
fprintf('Fix Tmin = %.1f degC, Tmax = %.1f degC at literature values.\n', Tmin, Tmax);
fprintf('Remaining parameters to estimate: A, C, M, a, b (p=5)\n');
%%
%[text] ## ratio of SSC pairs (all 7 params)
% Check for linear dependence between SSC columns
fprintf('\n--- SSC ratio analysis (all 7 params) ---\n');
for i = 1:p_all
    for j = i+1:p_all
        ratio_ij = Xp_all(:,i)./Xp_all(:,j);
        ratio_std = std(ratio_ij(isfinite(ratio_ij)));
        fprintf('  std(SSC_%s / SSC_%s) = %.6f', pnames_all{i}, pnames_all{j}, ratio_std);
        if ratio_std < 0.1
            fprintf('  <-- near-constant ratio (linearly dependent)\n');
        else
            fprintf('\n');
        end
    end
end
%%
%[text] ## ============ SSC for 5 estimable parameters (Tmin, Tmax fixed) ============
Tmin_fixed = Tmin;
Tmax_fixed = Tmax;

beta0(1) = A;
beta0(2) = C;
beta0(3) = M;
beta0(4) = a;
beta0(5) = b;
p = length(beta0);
pnames = {'A','C','M','a','b'};
%%
%[text] ## define the function that will be used for the forward problem (5 params)
fnameFOR = @(beta,t) gompertzFOR_5(beta, t, Tmin_fixed, Tmax_fixed);
%%
%[text] ## X' = SSC using initial guesses -- 5 estimable parameters
Xp = SSC_V3(beta0, xs, fnameFOR);
%%
%[text] ## plot X' for 5 parameters
clear h2
figure
hold on
set(gca, 'fontsize',18,'fontweight','bold');
for i = 1:p
    h2(i) = plot(xs(1:ns), Xp(1:ns,i), '-', 'color', cmap(i,:), 'LineWidth', 3);
end
ypred = fnameFOR(beta0, xs);
h2(p+1) = plot(xs, ypred, '--', 'color', cmap(p+1,:), 'LineWidth', 4);
legStr = cell(1, p+1);
for i = 1:p
    legStr{i} = [pnames{i}, '*\partialY/\partial(', pnames{i}, ')'];
end
legStr{p+1} = 'Y';
legend(h2, legStr, 'location', 'best');
xlabel('time (hr)');
ylabel('scaled sensitivity coefficient \beta_i*\partial(Y)/\partial(\beta_i), Y units');
title('SSC using initial guesses -- 5 estimable parameters');
grid on
%%
%[text] ## ratio of SSC pairs (5 params) -- check linear independence
fprintf('\n--- SSC ratio analysis (5 estimable params) ---\n');
for i = 1:p
    for j = i+1:p
        ratio_ij = Xp(:,i)./Xp(:,j);
        ratio_std = std(ratio_ij(isfinite(ratio_ij)));
        fprintf('  std(SSC_%s / SSC_%s) = %.6f', pnames{i}, pnames{j}, ratio_std);
        if ratio_std < 0.1
            fprintf('  <-- near-constant ratio (linearly dependent)\n');
        else
            fprintf('\n');
        end
    end
end
%%

%%
%[text] ## SSC summary table
fprintf('\n========== SSC Summary (5 estimable params, initial guesses) ==========\n');
fprintf('%-10s  %12s  %12s  %12s\n', 'Parameter', 'Max|SSC|', 'Min|SSC|', 'Mean|SSC|');
fprintf('%-10s  %12s  %12s  %12s\n', '---------', '--------', '--------', '---------');
for i = 1:p
    absSSC = abs(Xp(:,i));
    fprintf('%-10s  %12.6f  %12.6f  %12.6f\n', pnames{i}, max(absSSC), min(absSSC), mean(absSSC));
end
%%
%[text] ## Condition number of sensitivity matrix
condXp_all = cond(Xp_all);
condXp = cond(Xp);
fprintf('\n========== Condition Number ==========\n');
fprintf('  cond(X'') for 7 params: %.4f\n', condXp_all);
fprintf('  cond(X'') for 5 params: %.4f\n', condXp);
fprintf('  (should be < 1e6 for reliable estimation)\n');
%%
%[text] ## det(X''T * X'') -- estimability check
detXTX_all = det(Xp_all'*Xp_all);
detXTX = det(Xp'*Xp);
fprintf('\n========== det(X''T * X'') ==========\n');
fprintf('  7 params: %.6e\n', detXTX_all);
fprintf('  5 params: %.6e\n', detXTX);
fprintf('  (should be far from zero)\n');

%%
%[text] ## functions used in this code

%% Forward model -- ALL 7 parameters (for initial SSC analysis)
function logN = gompertzFOR_7(beta, t)
    global tTemp
    A=beta(1); C=beta(2); M=beta(3);
    a=beta(4); b=beta(5); Tmin=beta(6); Tmax=beta(7);
    tspan = t(:);
    T0 = interp1(tTemp(:,1), tTemp(:,2), tspan(1), 'linear', 'extrap');
    mu0 = secModel(T0, a, b, Tmin, Tmax);
    K0 = exp(-mu0*(tspan(1)-M));
    y0 = [A; K0];
    [~,Y] = ode45(@(t,y) gompODE(t,y,C,M,a,b,Tmin,Tmax), tspan, y0);
    logN = Y(:,1);
end

%% Forward model -- 5 parameters (Tmin, Tmax fixed)
function logN = gompertzFOR_5(beta, t, Tmin, Tmax)
    global tTemp
    A=beta(1); C=beta(2); M=beta(3);
    a=beta(4); b=beta(5);
    tspan = t(:);
    T0 = interp1(tTemp(:,1), tTemp(:,2), tspan(1), 'linear', 'extrap');
    mu0 = secModel(T0, a, b, Tmin, Tmax);
    K0 = exp(-mu0*(tspan(1)-M));
    y0 = [A; K0];
    [~,Y] = ode45(@(t,y) gompODE(t,y,C,M,a,b,Tmin,Tmax), tspan, y0);
    logN = Y(:,1);
end

%% ODE: d(logN)/dt = mu*K*C*exp(-K), dK/dt = -mu*K
function dy = gompODE(t, y, C, M, a, b, Tmin, Tmax)
    global tTemp
    K = y(2);
    T = interp1(tTemp(:,1), tTemp(:,2), t, 'linear', 'extrap');
    mu = secModel(T, a, b, Tmin, Tmax);
    dlogN = mu*K*C*exp(-K);
    dK    = -mu*K;
    dy = [dlogN; dK];
end

%% Secondary model: mu = a*(T-Tmin)^2*(1-exp(b*(T-Tmax)))
function mu = secModel(T, a, b, Tmin, Tmax)
    if T <= Tmin || T >= Tmax
        mu = 0;
    else
        mu = a*(T-Tmin)^2*(1-exp(b*(T-Tmax)));
    end
end

%[text] ## X' = scaled sensitivity coefficients using forward-difference
%[text] This is a forward problem with known approximate parameters
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
