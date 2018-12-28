%% Function to find the partition function for triatomic gases
function Qint = Q1(data,T)

% Data taken from
% https://physics.nist.gov/PhysRefData/Handbook/Tables/oxygentable5.htm
% Move around in the same website to find of all other elements

% Other constants 
global h c k

% Input from file
Te = data(:,1);
p  = data(:,2);

Qint = 0;
for n = 1:size(Te)
    
    Qel = p(n)*exp(-(Te(n))*h*c*100/(k*T));
    
    Qint = Qint + Qel;
    
end

end
