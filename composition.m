%% Program to find the equilibrium composition of air
% There are seven species in all namely N2 O2 NO NO+ N O e
close all; clear all; clc;

%% INPUT

P_tot = 90*1.01325e5;        % Total pressure
% T     = 2000;        % Temperature

XN2   = 0.79;       % Mole fraction N2
XO2   = 0.21;       % Mole fraction O2

P0_N2 = P_tot*XN2
P0_O2 = P_tot*XO2

data_N2  = dlmread('KpN2.dat');
data_O2  = dlmread('KpO2.dat');
data_NO  = dlmread('KpNO.dat');
data_NOe = dlmread('KpNO+.dat');

T_plot    = [];
Mole_frac = [];


    X0 = [P0_N2 P0_O2 0 0 0 0 0]';
%     X0 = [P0_N2 P0_O2 0 0 0]';

for T = 300:750:30000
    
    T_plot = [T_plot T];
    
    Kp_N2  = spline(data_N2(:,1),data_N2(:,2),T);
    Kp_O2  = spline(data_O2(:,1),data_O2(:,2),T);
    Kp_NO  = spline(data_NO(:,1),data_NO(:,2),T);
    Kp_NOe = spline(data_NOe(:,1),data_NOe(:,2),T);

    % Solve nonlinear equations
%     options = [];
    options = optimoptions('fsolve','TolFun',1e-9,'MaxFunEvals',15000,'Display','none');
    sol     = fsolve(@comp_rhs,X0,options,P0_N2,P0_O2,Kp_N2,Kp_O2,Kp_NO,Kp_NOe);
%     [sol,f]     = fsolve(@comp_rhs,X0,options,P0_N2,P0_O2,Kp_N2,Kp_O2,Kp_NO);
    
    % Mole fractions of each of the species
% sum(sol)
    
    sol(4) = abs(sol(4));
    sol(7) = abs(sol(7));
    X_n = sol/sum(sol);
    X0 = sol;
    Mole_frac = [Mole_frac X_n];
    
end

OUT = [T_plot' Mole_frac'];

%% Plotting

plot(T_plot,(Mole_frac),'-o')
grid on 
grid minor
title('Mole fraction at P = 90atm')
xlabel('Temperature')
ylabel('X_i')
legend('location','best','N2','O2','NO','NO+','N','O','e')

%% Save to file

% save comp_P_75atm.dat OUT -ASCII

