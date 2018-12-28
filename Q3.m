%% Function to find the partition function for triatomic gases
function Qint = Q3(data,T)

% Other constants 
global h c k

% Input from file
Te      = data(:,1);
p       = data(:,2);
A0      = data(:,4);
B0      = data(:,5);
C0      = data(:,6);
v1      = data(:,8);
v2      = data(:,9);
v3      = data(:,10);
Sigma   = data(:,11);
Sigma_r = data(:,12);

% Degeneracies in vibrational states
d1 = 1;
d2 = 1;
d3 = 1;

Qint = 0;
for n = 1:size(Te)

    Qrr = 0;
    Qho = 0;
    Qel = 0;
    
    % Electronic partition function
    Qel = p(n)*exp(-(Te(n))*h*c*100/(k*T));
    
    % Harmonic oscilator
    Qho = ((1-exp(-100*h*c*v1(n)/(k*T)))^(-d1))*...
           ((1-exp(-100*h*c*v2(n)/(k*T)))^(-d2))*...
           ((1-exp(-100*h*c*v3(n)/(k*T)))^(-d3));
    
    % Rigid rotor
    if Sigma_r(n) == 0
        Qrr = k*T/(100*Sigma(n)*h*c*B0(n));
    else
        Qrr = (pi/(A0(n)*B0(n)*C0(n)))^0.5 * (k*T/(100*h*c))^1.5;
    end
    
    Qint = Qint + Qel*Qrr*Qho;
    
end

end

    