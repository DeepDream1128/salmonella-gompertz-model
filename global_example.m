%% example of using global to share an array between script and function
global tTemp
tTemp=readmatrix('timeTemp.xls');
%first column of tTemp is time
%second column of tTemp is temperature
data=readmatrix('xyobs.xls');%these are the times and temperatures for the experimental data
xobs=data(:,1); %note that these times xobs are different from the times for temperature
yobs=data(:,2);

[beta,resids,J,COVB,mse] = nlinfit(xobs,yobs,fnameINV,beta0);

function y=model(beta,t)
global tTemp
[t,y]=ode45(@ff,tspan,beta(1)); %beta(1) is the initial value and is a parameter
    function dy=ff(t,y)
        T=interp1(tTemp(:,1),tTemp(:,2),t); %interpolates temperature at any time
        dy(1)=K*C*(T-35)*beta(2);
    end
end

