close all; clear all; clc;

global h c k
data  = xlsread('NO+');
T = 10000;
h     = 6.626068e-34;                      % Planck's constant - m^2 kg / s
c     = 299792458;                         % Speed of light - m/s
k     = 1.3806e-23;                        % Boltzmann constant - SI



Q = Q2approx(data,T)

Q = Q2(data,T)

cp = Cp(data,T)

% % 
% % %%
T_plot = [];
% Q = [];
cp = [];
for T = 100:3000:21100
    T
    T_plot = [T_plot T];
    cp = [cp Cp(data,T)];

end

% 

plot(T_plot,cp,'o')
grid on 
grid minor
title('Cp/R')














% 
% global h c
% 
% h  = 6.626068e-34;          % Units are m^2 kg / s
% Mu = 1.072988784712464e-26;            % kg
% c  = 299792458;           % m/s
% 
% %% Partition function of O2 - Spectroscopic data
% 
% % For first electronic state - T0 = 0
% We = 1580.19;
% WeXe = 11.98;
% WeYe = 0.04747;
% WeZe = -0.001273;
% D0 = 41260.08;
% g = 3;
% Be = 1.4378;
% De = 4.839e-6;
% ae = 1.593e-2;
% beta = 0;
% re = 1.2075e-8;
% 
% %% Pure vibration - V_max
% 
% E0 = 0.5*We - 0.25*WeXe + 0.125*WeYe + WeZe/16
% W0 = We - WeXe + 0.75*WeYe + 0.125*WeZe;
% W0X0 = WeXe - 1.5*WeYe - 1.5*WeZe;
% W0Y0 = WeYe + 2*WeZe;
% W0Z0 = WeZe;
% 
% p = [W0Z0 W0Y0 -W0X0 W0 -D0];
% r = roots(p);
% r = r(imag(r)==0);
% r = r(r>0);
% V_max = ceil(min(r));
% 
% %% Rotation states - J_max
% 
% V = 0
% Jguess = 100;
% Evib = E0 + W0*V - W0X0*V^2 + W0Y0*V^3 + W0Z0*V^4
% b = We*sqrt(2*pi^2*c*Mu/(100*De*h))
% Bv = Be - ae*(V+0.5)
% Dv = De - beta*(V+0.5)
% 
% % D0 = 42047;
%  x0 = [1e-3 Jguess];
%  options = optimoptions('fsolve','TolFun',1e-6,'MaxFunEvals',200,'Display','none');
%  sol = fsolve(@Jmax_eq,x0,options,D0,b,re,Mu,Evib,Bv,Dv);
%  J_max = sol(2)
%  
%  if J_max<0
%             J_max=0;
%         else
%             J_max=ceil(sol(2));
%         end
% 
% J_max

