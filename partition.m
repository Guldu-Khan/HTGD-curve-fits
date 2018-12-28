%% Program to save partition function data
close all; clear all; clc;

%% INPUT

global h c k
data  = xlsread('e');

h     = 6.626068e-34;                      % Planck's constant - m^2 kg / s
c     = 299792458;                         % Speed of light - m/s
k     = 1.3806e-23;                        % Boltzmann constant - SI

% m     = 4*data(1,13);
% m = 4.99186e-26;    % NO
% m = 2.33593e-26;        % N
% m = 2.6552e-26;         % O
m = 9.10938356e-31;         % e

 Qint = [];
 Q    = [];
 
 T_plot = [50 100 200 300 400 500 600 700 800 900 1000 2000 3000 5000 8000 10000 12500 15000 18000 20000 22500 25000 28000 30000]';

%  Q2(data,30000)
 
 for T = [50 100 200 300 400 500 600 700 800 900 1000 2000 3000 5000 8000 10000 12500 15000 18000 20000 22500 25000 28000 30000]
     
     T
     qint = Q1(data,T);
     qtr  = (2*pi*m*k*T/h^2)^1.5;
     Qint = [Qint qint];
     Q    = [Q qtr*qint];
     
 end
 
 OUT = [T_plot Qint' Q'];
 
 %% Plot 
 
 plot(T_plot,Qint,'ko');
 grid on
 grid minor
 title('Qint')
 
 %% Save data to file
 
%  save Qe.dat OUT -ASCII
    
    