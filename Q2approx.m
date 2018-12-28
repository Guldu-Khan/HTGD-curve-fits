%% Program for partition function for N2 using Drellishak approx
function Qint = Q2approx(data,T)

% Input from xlsx file

% All spectroscopic data required is first input using the required file.
% All the data below are spectroscopic constants.

Te    = data(:,1);
g     = data(:,2);
We    = data(:,3);
WeXe  = data(:,4);
WeYe  = data(:,5);
WeZe  = data(:,6);
D0    = data(:,7);
Be    = data(:,8);
ae    = data(:,9)*1e-2;
De    = data(:,10)*1e-6;
beta  = data(:,11)*1e-6;
re    = data(:,12)*1e-8;
Mu    = data(1,13);             % kg 
sigma = data(1,14);             % Homonuclear/Heteronuclear

% Other constants 
global h c k

%% Calculating the internal partition function.

% N_max = size(Te,1);     % N2
N_max = 4;
Qint = 0;

for n = 1:N_max
    
    Qel = g(n)*exp(-(Te(n))*h*c*100/(k*T));
    
    Bel  = Be(n);
    dele = ae(n)/Bel;
    sig  = 100*(1-0.5*dele)*Bel*h*c/(k*T);
    xe   = WeXe(n)/We(n);
    u    = 100*(1-2*xe)*We(n)*h*c/(k*T);
    gam  = Bel/We(n);
    big  = (1/(sig*(1-exp(-u))))*(1+(sig/3)+(sig^2/15)+(8*gam*gam/sig)-(dele/(1-exp(u)))+(2*xe*u/(1-exp(u))^2));
    Qint = Qint + (1/sigma)*Qel*big;
    
end

end




