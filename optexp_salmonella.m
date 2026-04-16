%% optimal exptl design
% this is a forward problem only reproduce Beck and Arnold, p. 442, Fig. 8.10 
% Salmonella Gompertz model: 5 estimated parameters [A,C,M,a,Tmax]
% Secondary model: mu = a*(T-Tmin)^2*(1-exp(b*(T-Tmax)))
%% set up times

close all
clear
global tTemp
s = settings;
s.matlab.fonts.editor.code.Size.TemporaryValue = '14pt';

% --- Load temperature data (needed by the Gompertz ODE model) ---
tempData = readmatrix('Salmonella sin growth Temps.xlsx');
tTemp_hr = tempData(:,1) / 60;   % minutes -> hours
Tobs     = tempData(:,2);
tTemp    = [tTemp_hr, Tobs];

t=linspace(0, max(tTemp_hr), 300)';%time grid spanning temperature data
m=length(t);
%% parameter values -- OLS estimates for [A, C, M, a, Tmax]

beta(1)=2.3581;   % A  (initial log10 CFU/mL)
beta(2)=11.0623;  % C  (asymptotic growth range)
beta(3)=5.5406;   % M  (inflection time, hr)
beta(4)=0.0003;   % a  (secondary model coefficient)
beta(5)=45.7230;  % Tmax (upper temperature limit, degC)

Tmin_fixed = 6.0;       % fixed from literature
b_fixed    = 0.275;      % fixed from literature
pnames = {'A','C','M','a','T_{max}'};
%% Call and plot the function

fname=@(b,t) gompertzFOR(b,t,Tmin_fixed,b_fixed);
Y=fname(beta,t);
plot(t,Y)
xlabel('time (hr)'); ylabel('log_{10}N'); title('Gompertz model prediction');
%Make sure Y and t have exactly the same dimensions
%Do not allow Y to be a row and t to be a column.  Make both of them
%columns
%% Compute C11, C12,...delta for opt exptl design

%Xp are the scaled sensitivity coefficients
Xp=SSC_V4(beta,t,fname);%compute scaled sensitivity coefficients as a cell array
p=length(beta);
%% compute entire C matrix

for i=1:p 
    for j=1:p
        intgrnd=Xp{i}.*Xp{j}; %integrand for Eq. 8.3.5, Beck and Arnold, p. 434
        C{i,j}=(1./t).*cumtrapz(t,intgrnd);  % trapezoidal rule integral of Eq. 8.3.5, Beck and Arnold, p. 434
        clear intgrnd
    end
end
%% to compute delta, must set up the C matrix for each time

%extract C into a 3-D matrix we call "CC"
for i=1:p
    for j=1:p
        CC(i,j,:)=C{i,j}; %gives a 3D matrix that is m (depth) in time
    end
end
%% For your project, you must set all CC(:,:,1) (first page) values to the appropriate
% number, based on your project

CC(:,:,1)=0;%beginning values at time = 0 are zero, except if the parameter is Y(0) = initial values
CC(1,1,1)=1;%except for C11, for the initial value (A = logN(0))
delta(1)=0; %determinant at time zero = 0.
for k=2:m
    delta(k)=det(CC(:,:,k));%delta at each time
end
% delta=sqrt(delta);%converts units of delta to same as units for C
%% convert 3D to 2D for plotting

%C matrix is symmetrical, so need only one half
% "Cp" is "C for plotting"
% Cp(i) holds Cii for the i-th parameter
for i=1:p
    Cp(i,:)=CC(i,i,:);
end
%% plot the results

figure
hold on
fac=32;%factor for plotting
cmap = ['r' 'g' 'b' 'c' 'k']';
for i=1:p
    h(i)=plot(t,Cp(i,:),'-','color',cmap(i,:),'linewidth',3);
end
f_delta = round(max(max(Cp(:,2:end))) / max(delta(2:end)));
if f_delta < 1, f_delta = 1; end
h(p+1)=plot(t,f_delta*delta,'-^','color',[0.5 0 0.5],'linewidth',3,'MarkerSize',5,'MarkerIndices',1:20:m);
set(gca, 'fontsize',28,'fontweight','bold');
xlabel('time (hr)')
legstr = pnames;
legstr{end+1} = sprintf('%d\\Delta', f_delta);
legend(legstr,'Location','Best')
grid on
%% functions

function logN = gompertzFOR(beta, t, Tmin, b)
    global tTemp
    A = beta(1); C = beta(2); M = beta(3);
    a = beta(4); Tmax = beta(5);
    tAll = t(:);
    [tUniq, ~, ic] = unique(tAll);
    y0 = A;
    [~,Y] = ode45(@(tt,y) gompODE(tt,y,C,M,a,b,Tmin,Tmax), tUniq, y0);
    logN = Y(ic,1);
end

function dy = gompODE(t, y, C, M, a, b, Tmin, Tmax)
    global tTemp
    T = interp1(tTemp(:,1), tTemp(:,2), t, 'linear', 'extrap');
    mu = a*(T - Tmin)^2 * (1 - exp(b*(T - Tmax)));
    K = exp(-mu*(t - M));
    dy = mu*K*C*exp(-K);
end

function Xp=SSC_V4(beta,x,yfunc)
%Computes scaled sensitivity coefficients =Xp, nxp matrix
%Returns scaled sensitivity matrix in a cell array
%can have k dependent variables that are stacked in a column vector
%all y1s, then all y2s, ...last are yks
%n is the number of data
%p is the number of parameters
%Xp1 = dY/dbeta1~[y(beta1(1+d), beta2,...betap) - y(beta1,
%beta2,...betap)]/d...
%d is the arbitrary delta
%beta is the p x 1 parameter vector
%yhat is nx1 vector, the y values when only one parameter has been successively perturbed by d
%ypred is nx1 vector,  the y values when parameters are set to beta
%betain is px1 vector, the parameter values with only one parameter perturbed by d
%x are the independent variables (can be one or more independent variables)
%yfunc is a function (m file or an anonymous) defined by the user outside
%of this file
%% X' = scaled sensitivity coefficients using forward-difference
% This is a forward problem with known approximate parameters

d=0.001;
ypred=yfunc(beta,x);
for i = 1:length(beta)  %scaled sens coeff for forward problem
    betain = beta; %reset beta
    betain(i) = beta(i)*(1+d);
    yhat{i} = yfunc(betain,x);
    Xp{i} = (yhat{i}-ypred)/d;%scaled sens coeff for ith parameter
end
end
