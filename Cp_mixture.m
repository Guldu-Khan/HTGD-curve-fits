%% Program to find Cp of air.
close all; clear all; clc;

%% INPUT

% Gas constant
R = 8.314;

% Composition
comp = dlmread('comp_P_90atm.dat');

% Cp individual
data_N2  = dlmread('CpN2.dat');
data_O2  = dlmread('CpO2.dat');
data_NO  = dlmread('CpNO.dat');
data_NO1 = dlmread('CpNO+.dat');
data_N   = dlmread('CpN.dat');
data_O   = dlmread('CpO.dat');
data_e   = dlmread('Cpe.dat');

%% Plot each one

T_p = comp(:,1);
n_p = comp(:,2:end);   

T_plot = 200:100:30000;
np_1 = spline(T_p,n_p(:,1),200:100:30000);
np_2 = spline(T_p,n_p(:,2),200:100:30000);
np_3 = spline(T_p,n_p(:,3),200:100:30000);
np_4 = spline(T_p,n_p(:,4),200:100:30000);
np_5 = spline(T_p,n_p(:,5),200:100:30000);
np_6 = spline(T_p,n_p(:,6),200:100:30000);
np_7 = spline(T_p,n_p(:,7),200:100:30000);

N_p = [np_1' np_2' np_3' np_4' np_5' np_6' np_7'];
figure(1)
plot(T_plot,N_p);
grid on;
grid minor;
title('Composition')
xlabel('Temperature')
ylabel('Mole fraction')
legend('N2','O2','NO','NO+','N','O','e')

%% Find Cp mixture

T_out  = [];
Cp_mix = [];

for T = 300:750:30000
    
    T_out  = [T_out T];
    
    % Mole fraction
    nN2  = spline(comp(:,1),comp(:,2),T);
    nO2  = spline(comp(:,1),comp(:,3),T);
    nNO  = spline(comp(:,1),comp(:,4),T);
    nNOe = spline(comp(:,1),comp(:,5),T);
    nN   = spline(comp(:,1),comp(:,6),T);
    nO   = spline(comp(:,1),comp(:,7),T);
    
    % Mass fraction
    M    = 28*nN2 + 32*nO2 + 30*nNO + 30*nNOe + 14*nN + 16*nO;
    
    CpN2 = spline(data_N2(:,1),data_N2(:,2),T)*R/28;
    CpO2 = spline(data_O2(:,1),data_O2(:,2),T)*R/32;
    CpNO = spline(data_NO(:,1),data_NO(:,2),T)*R/30;
    CpNOe= spline(data_NO1(:,1),data_NO1(:,2),T)*R/30;
    CpN  = spline(data_N(:,1),data_N(:,2),T)*R/14;
    CpO  = spline(data_O(:,1),data_O(:,2),T)*R/16;
    
    Cp = (CpN2*28*nN2 + CpO2*32*nO2 + CpNO*30*nNO + CpNOe*30*nNOe + CpN*14*nN + CpO*16*nO)/M;
    
    Cp_mix = [Cp_mix Cp];
    
end

OUT = [T_out' Cp_mix'];

%% Plot Cp of air

figure(2)
plot(T_out,Cp_mix,'-ko')
grid on
grid minor
title('Cp mixture at P=90atm')
xlabel('Temperature')
ylabel('Cp')

%% Save data to file

%  save Cpair_90atm.dat OUT -ASCII

    
    
    
