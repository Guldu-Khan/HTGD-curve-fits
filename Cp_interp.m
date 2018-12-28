%% Program to use Q data and form Cp
close all; clear all; clc;

%% INPUT

data = dlmread('QO2.dat');

T    = data(:,1);
Qint = data(:,2);
% Q    = data(:,3);

Cp = [];
T_plot = [];
dT = 0.25;
Q_plot = [];

for j = 400:200:30000
    
    T_plot = [T_plot j];
    Q_plot = [Q_plot spline(T,Qint,j)];
%     
    lnQT2   = log(spline(T,Qint,j+dT));
    lnQT0   = log(spline(T,Qint,j-dT));
    lnQT1   = log(spline(T,Qint,j));
    dlnqdT  = (lnQT2 - lnQT0) /(2*dT);
    d2lnqdT = (lnQT2+lnQT0-2*lnQT1)/dT^2;
%     
%     lnQT2   = log(spline(T,Q,j+dT));
%     lnQT0   = log(spline(T,Q,j-dT));
%     lnQT1   = log(spline(T,Q,j));
%     dlnqdT  = (lnQT2 - lnQT0) /(2*dT);
%     d2lnqdT = (lnQT2+lnQT0-2*lnQT1)/dT^2;
    
    Cpint   = 2*j*dlnqdT + j^2*d2lnqdT;

    cp      = Cpint + 2.5;
    
    Cp = [Cp cp];
    
end

OUT = [T_plot' Cp'];

%% Data from Capitelli for plotting

% Tcap = [200 400 600 800 1000 2000 4000 6000 8000 10000 15000 17000 20000 25000 30000];
% Ccap = [1.001 1.018 1.121 1.281 1.433 1.827 2.017 2.106 2.312 2.925 5.438 5.617 4.858 3.040 1.841];
% Qcap = [3.512e1 7.011e1 1.055e2 1.423e2 1.816e2 4.311e2 1.259e3 ...
%     2.546e3 4.33e3 6.701e3 1.729e4 2.498e4 4.261e4 9.345e4 1.754e5];  % N2  
Tcap = [200 400 600 800 1000 2000 4000 6000 8000 10000 15000 20000 25000 30000];
Ccap = [1.003 1.121 1.360 1.558 1.695 2.046 2.518 2.792 2.813 2.520 1.554 .9960 .7096 .5437];
Qcap = [1.464e2 2.937e2 4.498e2 6.239e2 8.209e2 2.207e3 7.452e3 1.719e4  ...
        3.302e4 5.624e4 1.474e5 2.743e5 4.195e5 5.706e5]; % O2

%% Plot Cp

figure(1)
plot(T_plot,Cp,'k',Tcap(1:end-3),Ccap(1:end-3)+2.5,'ro');
grid on;
grid minor;
title('Specific heat of O2')
xlabel('Temperature')
ylabel('Cp/R')
legend('location','best','Present computation','Capitelli')

figure(2)
plot(T_plot,Q_plot,'k',Tcap,Qcap,'ro')
grid on
grid minor
title('Internal partition function of O2')
xlabel('Temperature')
ylabel('Q_{int}')
legend('location','best','Present computation','Capitelli')

%% Save data to file

%  save Cpe.dat OUT -ASCII
    