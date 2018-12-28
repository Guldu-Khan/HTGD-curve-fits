%% Program to estimate parameters for curve fitting
close all; clear all; clc;

%% INPUT

global h c k data

h     = 6.626068e-34;                      % Planck's constant - m^2 kg / s
c     = 299792458;                         % Speed of light - m/s
k     = 1.3806e-23;                        % Boltzmann constant - SI

data   = dlmread('Cpair_1atm.dat');

%% Plot Cp

T_p    = data(:,1);
Cp_p   = data(:,2);

T_plot = 200:1:30000;
Cp_1 = spline(T_p,Cp_p,200:1:30000);

plot(T_plot,Cp_1,'-k');
grid on;
grid minor;
title('Cp/R')
xlabel('T')
ylabel('Cp')

%% Parameter estimation

% control_pts = [100 2500 8000 18000 22000 30000]';   % For N2
% control_pts = [100 1500 3000 5000 8000 15000 30000]'; % For O2
control_pts = [200 1000 3500 5000 6300 8500 10000 12500 20000 30000]'; % For air 1 atm

parameter   = [];

for j = 1:(size(control_pts,1)-1)
    
    T = linspace(control_pts(j),control_pts(j+1),10);
    X0 = zeros(7,1);
    ObjectiveFunction = @curvefit_cost;
    options = optimset('MaxIter',25000,'MaxFunEvals',25000,'TolFun',1e-9,'TolX',1e-9);

    [X, fval, exitFlag,output] = fminsearch(ObjectiveFunction,X0,options,T);

    parameter = [parameter X];
    
end

ctr = control_pts';

%% Save parameter values to file

% save par_air_1atm.dat ctr parameter -ASCII

%%
% T = 100;
% Cp = spline(data(:,1),data(:,2),T)
% Cp_fit = X(1)/T^2 + X(2)/T + X(3) + X(4)*T + X(5)*T^2 + X(6)*T^3 + X(7)*T^4