%% Program to find Kp of reactions
close all; clear all; clc;

%% INPUT

% D0 = 1.6022e-19*9.753940;       % Input for N2 dissociation
D0 = 1.6022e-19*5.116420;          % O2
% D0 = 1.6022e-19*6.506;          % NO
% D0 = 1.6022e-19*10.8506;          % NO+ --> N + O+
% D0 = 1.6022e-19*9.2642;             % NO --> NO+ + e-

% m = 2.335934e-26/2;               % N2 % ma*mb/mab;
m = 1.32815e-26;                    % O2
% m = 1.24281E-26;                    % NO / NO+
% m = 9.10938356e-31;                 % e

global h c k

h     = 6.626068e-34;                      % Planck's constant - m^2 kg / s
c     = 299792458;                         % Speed of light - m/s
k     = 1.3806e-23;                        % Boltzmann constant - SI

data_A = dlmread('QO.dat');
% data_B = dlmread('QO.dat');
data_AB = dlmread('QO2.dat');


T    = data_A(:,1);
QA   = data_A(:,2);
% QB   = data_B(:,2);
QAB   = data_AB(:,2);

Kp     = [];
T_plot = [];

for j = 1000:100:30000
    
    T_plot = [T_plot j];
    
    qa  = spline(T,QA,j);
%     qb  = spline(T,QB,j);
    qab = spline(T,QAB,j);
    
    kp  = exp(-D0/(k*j))*(qa^2/(qab))*(2*pi*m*k*j/h^2)^1.5*k*j;
%     kp  = exp(-D0/(k*j))*(qa*qb/(qab))*(2*pi*m*k*j/h^2)^1.5*k*j;
        
    Kp = [Kp kp];
    
end

OUT = [T_plot' Kp' log10(Kp')];

%% Bartlem data

Tbar = [1000 2000 4000 6000 8000 10000];
% lKbar = [-3.8e1 -1.3077e1 -4.96e-1 3.7552 5.93537 7.28996]; % N2
lKbar = [-1.46132e1 -1.35267 5.3431 7.58653 8.70884 9.3769]; % O2

%% Plotting

plot(T_plot,log10(Kp),'k',Tbar,lKbar,'ro');
grid on;
grid minor;
title('O2 --> 2O')
xlabel('Temperature')
ylabel('log Kp')
legend('location','best','Present computation','Bartlem');

%% Save to file

% save KpN2.dat OUT -ASCII
    
    
    
    

