%% example of nlinfit and ode45 for two equations using std dev as weighting factors

%Chapra 22.19 Inverse Problem on the following equations:
%k=exp(-beta(2)*(1/(T+273)-1/(95+273.15)))   k is a rate constant
%dC/dt = -kC    C is the first dependent variable
%dT/dt = beta(1)*k*C - beta(3)(T-20)   T is the second dependent variable
%% Preprocessing

s = settings;
s.matlab.fonts.editor.code.Size.TemporaryValue = '14pt';
clear  %clear all variables
close all
format compact
data =readmatrix('Seborg20_19_data2.xlsx');
%% initial guesses
% 

beta0(1)=8e2;
beta0(2)=1e3;
beta0(3)=10; 
beta0(4)=1.5; %C(0)
beta0(5)=15;  %T(0)
p=length(beta0);
%% Experimental data

xobs=data(:,1);x2=xobs; %x2 for plotting
n=length(xobs);
yobsC=data(:,2);
yobsT=data(:,3);
stdC=.07;stdT=3; %estimates of the standard deviation for each dependent variable 
% stdC=1;stdT=1;% if we do not use weighted least squares
yobsC2=yobsC/stdC; yobsT2=yobsT/stdT; %regress on these scaled yobs (weighted least squares)
yobs=[yobsC2; yobsT2];%yobs for weighted least squares
%% define the function that will be used for the forward problem

fnameFOR=@Ex22_19_FOR;
%% plot ypred with initial parameter guesses on data

%to make sure guesses are reasonable
xs=linspace(min(xobs),max(xobs),100)'; %xs are the times for SSCs to make a smooth curve.
ns=length(xs);%length of xs for plotting
ypredInit=fnameFOR(beta0,xs);
ypredCInit=ypredInit(1:ns);
ypredTInit=ypredInit(ns+1:2*ns);
%% Output figure for concentration

%make vector for plotting x-axis
xp=xs(:,1);
figure
set(gca, 'fontsize',14,'fontweight','bold');
plot(xp, ypredCInit, '-','linewidth',2.5)
hold on %allows plotting more than one line on chart
plot(xobs, yobsC, 's','markersize',8,'markerfacecolor','r') 
xlabel('time (min)','fontsize',14,'fontweight','bold'); 
ylabel('Concentration gmol/L','fontsize',14,'fontweight','bold')
%% Output figure for temperature

figure
set(gca, 'fontsize',14,'fontweight','bold');
plot(xp, ypredTInit, '-','linewidth',2.5)
hold on
plot(xobs, yobsT, 's','markersize',8,'markerfacecolor','b'); 
xlabel('time (min)','fontsize',14,'fontweight','bold'); 
ylabel('Concentration gmol/L','fontsize',14,'fontweight','bold')
%% Scaled sensitivity coefficients before running inverse problem

%to determine a) which parameters can be estimated, and b) which parameters
%wil be most accurate (lowest relative error)
xs=linspace(min(xobs),max(xobs),100)'; %xs are the times for SSCs to make a smooth curve.
%make xs a column.
ns=length(xs);%length of xs for plotting
fnameFOR=@Ex22_19_FOR;
Xp=SSC_V3(beta0,xs,fnameFOR);
%% plot X' for C

%plot for C
cmap = ['r' 'g' 'b' 'c' 'k'  'm' 'y' ]';
figure
hold on
set(gca, 'fontsize',14,'fontweight','bold');
%plot C vs t to know the total span
ypred=fnameFOR(beta0,xs);%this computes ypred for both C (1:100) and T (101:200)
h2(1)=plot(xs(1:ns),ypred(1:ns),'-','color',cmap(1,:),'LineWidth',2); %plot the predicted C to compare to SSCs
for i=1:p
    h2(i+1) = plot(xs(1:ns),Xp(1:ns,i),'-','color',cmap(i+1,:),'LineWidth',2);
end
legend('C','\beta_1*\partialC/\partial\beta_1','\beta_2*\partialC/\partial\beta_2',...
    '\beta_3*\partialC/\partial\beta_3','C_o*\partialC/\partialC_o'...
    ,'T_o*\partialC/\partialT_o')
xlabel('time'); ylabel('scaled sensitivity coefficient for C, gmol/L');
grid on
%can check correlation between any 2 betas by dividing and plotting them
figure
rat1=Xp(1:ns,1)./Xp(1:ns,3);
figure
plot(xs(1:ns),rat1)
axis([0 5 -20 20]);
grid on
%% plot X' for T

figure
hold on
set(gca, 'fontsize',14,'fontweight','bold');
%plot T vs t to know the span
h3(1)=plot(xs(1:ns),ypred(ns+1:2*ns),'-','color',cmap(1,:),'LineWidth',2); %plot the predicted T to compare to SSCs
for i=1:p
    if i==1
        lw =3; %to make SSC1 thicker because it's on top of SSC4
    else
        lw=2;
    end
    h3(i+1) = plot(xs(1:ns),Xp(ns+1:2*ns,i),'-','color',cmap(i+1,:),'LineWidth',lw);
end
xlabel('time'); ylabel('scaled sensitivity coefficient for T, ^oC');
legend('T','\beta_1*\partialC/\partial\beta_1','\beta_2*\partialC/\partial\beta_2',...
    '\beta_3*\partialC/\partial\beta_3','C_o*\partialC/\partialC_o'...
    ,'T_o*\partialC/\partialT_o')
grid on
%can check correlation between any 2 betas by dividing and plotting their SSCs
figure
rat=Xp(101:200,1)./Xp(101:200,3);
figure
plot(xs(1:ns),rat)
axis([0 5 -10 5]);
grid on
%% nlinfit returns parameters, residuals, Jacobian (sensitivity coefficient matrix),

%covariance matrix, and mean square error.  ode45 is solved many times
%iteratively
fnameINV=@Ex22_19_INV;
xobs(1,2)=stdC;xobs(1,3)=stdT;%send the y stdev into the function for regression
[beta,resids,J,COVB,mse] = nlinfit(xobs, yobs,fnameINV, beta0);
rmse=sqrt(mse) %mean square error = SS/(n-p) total for weighted least squares
n=length(xobs);     nn=n(1); p=length(beta); 
beta
condX=cond(J) %must be < 1 million
detXTX=det(J'*J) % must not be near zero, the larger, the better
%rmse for each scaled dependent variable
% rC=resids(1:n); rT=resids(n+1:2*n);
% rmseC=sqrt(rC'*rC/(n-1))
% rmseT=sqrt(rT'*rT/(n-1))
%% R is the correlation matrix for the parameters, sigma is the standard error vector

[R,sigma]=corrcov(COVB)
relerr=sigma'./beta
%% confidence intervals for parameters

ci=nlparci(beta,resids,J)
%% computed ypredicted

ypred=Ex22_19_FOR(beta,xs);
ypredC=ypred(1:ns);
ypredT=ypred(ns+1:2*ns);
%% mean of the residuals

meanr=mean(resids)
%% Output figure for concentration

xobs=x2; %make xobs back into a 1-column vector for plotting
xp=xs(:,1); %create a plotting time vector
figure 
hold on
set(gca, 'fontsize',14,'fontweight','bold');
plot(xp, ypredC, '-','linewidth',2.5)
plot(xobs, yobsC, 's','markersize',8) 
xlabel('time (min)','fontsize',14,'fontweight','bold'); 
ylabel('Concentration gmol/L','fontsize',14,'fontweight','bold')
grid on
%% Figure for Temperature

figure
hold on
set(gca, 'fontsize',14,'fontweight','bold');
plot(xp, ypredT, '-r','linewidth',2.5)
plot(xobs, yobsT, 'or','markersize',8); %predicted T values
xlabel('time (min)','fontsize',14,'fontweight','bold'); ylabel('Temperature ^oC','fontsize',14,'fontweight','bold')
set(gca, 'fontsize',14,'fontweight','bold');
grid on
%% residual scatter plot

x3=[xobs; xobs;];
figure
hold on
h4(1)=plot(x3(1:n), resids(1:n), 'square','Markerfacecolor', 'b','markersize',10);
h4(2)=plot(x3(n+1:2*n), resids(n+1:2*n), 'o','Markerfacecolor', 'r');
YLine = [0 0]; 
XLine = [0 max(xobs)];
plot (XLine, YLine,'R'); %plot a straight red line at zero
ylabel('Observed y/\sigma - Predicted y/\sigma','fontsize',14,'fontweight','bold')
xlabel('time (min)','fontsize',14,'fontweight','bold')
legend(h4,'C','T')
%% number of runs = number of times moving from one residual to the next crosses zero

% rescross=resids(2:n).*resids(1:n-1);%multiply each pair of residuals
% res_sign=sign(rescross);%get the sign of each multiplied pair
% count=0;
% for i=1:n-1
%     if res_sign(i)<0 %if product of pair is < 0, that's a run
%         count=count+1;
%     end
% end
% fprintf('number of runs = %5.2f\n',count);
% minrun=(n+1)/2; %count should be >=minrun
% fprintf('Minimum required number of runs = %5.2f\n',minrun);
% This algorithm handles replicates and no replicates
%set up matrix with col1 = x, col2 = resids
x3=[xobs;xobs];%put xobs into x3
xResids = [x3 resids];
% sort residuals on time
xResidsSort=sortrows(xResids); %sorted xResids
x=xResidsSort(:,1); %sorted x
residsSort=xResidsSort(:,2); %sorted resids

count=0; %number of runs
countRep=1; %tracks number of reps 
countNeg=0; countPos=0;
resSign(countRep)=sign(residsSort(1)); %sign of first residual
for i=2:2*ns(end) %nS(end) is the total number of resids
    if x(i)==x(i-1) %if there is a replicate
        countRep=countRep+1;%increase counter for replicates
        resSign(countRep)=sign(residsSort(i));
    else %no replicate
        %count number of positive and negative residuals
        for j=1:countRep
            if resSign(j) < 0
                countNeg=countNeg+1;
            elseif resSign(j) > 0
                countPos=countPos+1;
            end
        end
        count=count+min(countNeg, countPos);
        %Check sign of product of current resid and previous resid
        rescross=residsSort(i)*residsSort(i-1); %product of residuals
        resSign(j+1)=sign(rescross);
        if resSign(j+1) < 0
            count=count+1;
        end
        clear resSign
        countRep=1;%reset
        resSign(countRep)=sign(residsSort(i));%set resSign(1)=sign of current resids
        countNeg=0; countPos=0;
    end
end
fprintf('number of runs = %5.2f\n',count);
minrun=(ns(end)+1)/2; %count should be >=minrun
fprintf('Minimum required number of runs = %5.2f\n',minrun);
%% residual histogram

figure
h=histogram(resids);
hold on
set(gca, 'fontsize',14,'fontweight','bold');
xlabel('Observed y/\sigma - Predicted y/\sigma','fontsize',16,'fontweight','bold')
ylabel('Frequency','fontsize',16,'fontweight','bold')
%% Scaled sensitivity coefficients

Xp=SSC_V3(beta,xs,fnameFOR);
%% plot X' for each dependent variable

%plot for C
cmap = ['r' 'g' 'b' 'c' 'k'  'm' 'y' ]';
figure
hold on
set(gca, 'fontsize',14,'fontweight','bold');
%plot C vs t to know the total span
ypred=fnameFOR(beta,xs);
h2(1)=plot(xs(1:ns),ypred(1:ns),'-','color',cmap(1,:),'LineWidth',2); %plot the predicted C to compare to SSCs
for i=1:p
    h2(i+1) = plot(xs(1:ns),Xp(1:ns,i),'-','color',cmap(i+1,:),'LineWidth',2);
end
legend('C','\beta_1*\partialC/\partial\beta_1','\beta_2*\partialC/\partial\beta_2',...
    '\beta_3*\partialC/\partial\beta_3','C_o*\partialC/\partialC_o'...
    ,'T_o*\partialC/\partialT_o')
xlabel('time'); ylabel('scaled sensitivity coefficient for C, gmol/L');
grid on
%% plot X' for T

figure
hold on
set(gca, 'fontsize',14,'fontweight','bold');
%plot T vs t to know the span
h3(1)=plot(xs(1:ns),ypred(ns+1:2*ns),'-','color',cmap(1,:),'LineWidth',2); %plot the prediced T to compare to SSCs
for i=1:p
    h3(i+1) = plot(xs(1:ns),Xp(ns+1:2*ns,i),'-','color',cmap(i+1,:),'LineWidth',2);
end
xlabel('time'); ylabel('scaled sensitivity coefficient for T, ^oC');
legend('T','\beta_1*\partialC/\partial\beta_1','\beta_2*\partialC/\partial\beta_2',...
    '\beta_3*\partialC/\partial\beta_3','C_o*\partialC/\partialC_o'...
    ,'T_o*\partialC/\partialT_o')
grid on
%% functions

function Xp=SSC_V3(beta,x,yfunc)
%Computes scaled sensitivity coefficients =Xp, nxp matrix
%can have k dependent variables that are stacked in a column vector
%all y1s, then all y2s, ...last are yks
%n is the number of data
%p is the number of parameters
%Xp1 = dY/dbeta1~[y(beta1(1+d), beta2,...betap) - y(beta1,
%beta2,...betap)]/d...
%d is the arbitrary delta, usually 0.001
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
    SSC{i} = (yhat{i}-ypred)/d;%scaled sens coeff for ith parameter
    Xp(:,i)=SSC{i}; %extract from cell array to 2D array
end
end

function y = Ex22_19_FOR(beta,t)
%function used for forward problem
%Used for SSC and plotting
%t column 1 are the times
%k=exp(-beta(2)*(1/(T+273.15)-1/(95+273.15))))
%y1 is C
%y2 is T
% global y0 %seqsort;
tspan=t(:,1); %we want y at every t
[t,y]=ode45(@batch,tspan,[beta(4) beta(5)]);

    function dy = batch(t,y) %function that computes the dydt
        %Chapra 20.19 solve the 2 ODEs
        k=@(T)exp(-beta(2)*(1/(T+273.15)-1/(95+273.15)));
        dy(1)=-k(y(2))*y(1);
        dy(2)=beta(1)*k(y(2))*y(1)-beta(3)*(y(2)-20);
        dy=dy';%places y values into a column
    end

% after the ode45, rearrange the n-by-2 y matrix into a 2n-by-1 matrix 
y1=y(:,1); %C=concentration, predicted values
y2=y(:,2); %T=temperature, predicted values
y=[y1;y2];%put the y's into a column
end

function y = Ex22_19_INV(beta,t)
%example of a function within a function to access and estimate parameters
%from Seborg 2.17
%t column 1 are the times
%t(1,2)=stdC  t(1,3)=stdT
%nlinfit will call this function from call_Ex22_19.m
%predicted ys are divided by the stdev for each dep variabel (weighted
%least squares = WLS)
%k=exp(-beta(2)*(1/(T+273.15)-1/(95+273.15))))
%y1 is C
%y2 is T
% global y0 %seqsort;
tspan=t(:,1); %we want y at every t
stdC=t(1,2);stdT=t(1,3);%stdev for each dep var
[t,y]=ode45(@batch,tspan,[beta(4) beta(5)]);

    function dy = batch(t,y) %function that computes the dydt
        %Chapra 20.19 solve the 2 ODEs
        k=@(T)exp(-beta(2)*(1/(T+273.15)-1/(95+273.15)));
        dy(1)=-k(y(2))*y(1);
        dy(2)=beta(1)*k(y(2))*y(1)-beta(3)*(y(2)-20);
        dy=dy';
    end

% after the ode45, rearrange the n-by-2 y matrix into a 2n-by-1 matrix
% and send that back to nlinfit. 
y1=y(:,1)/stdC; y2=y(:,2)/stdT;%predicted values
y=[y1;y2];%put the y's into a column
end
%% Comments: beta4 is the largest SSC and most uncorrelated, so it has lowest rel error
% beta5 has the larges rel error because SSC is so small

%the larger the SSC (if uncorrelated), the smaller the rel error