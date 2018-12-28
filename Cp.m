% Function to compute the Cp of species using the Partition function 
function cp = Cp(data,T)

dT = 0.25;

% For diatomic approximate
% lnQT2 = log(Q2approx(data,T+dT));
% lnQT0 = log(Q2approx(data,T-dT));
% lnQT1 = log(Q2approx(data,T));
% dlnqdT  = (lnQT2 - lnQT0) /(2*dT);
% d2lnqdT = (lnQT2+lnQT0-2*lnQT1)/dT^2;

% For diatomic
lnQT2 = log(Q2(data,T+dT));
lnQT0 = log(Q2(data,T-dT));
lnQT1 = log(Q2(data,T));
dlnqdT  = (lnQT2 - lnQT0) /(2*dT);
d2lnqdT = (lnQT2+lnQT0-2*lnQT1)/dT^2;


% For triatomic
% lnQT2 = log(Q3(data,T+dT));
% lnQT0 = log(Q3(data,T-dT));
% lnQT1 = log(Q3(data,T));
% dlnqdT  = (lnQT2 - lnQT0) /(2*dT);
% d2lnqdT = (lnQT2+lnQT0-2*lnQT1)/dT^2;

% For monoatomic
% lnQT2 = log(Q1(data,T+dT));
% lnQT0 = log(Q1(data,T-dT));
% lnQT1 = log(Q1(data,T));
% dlnqdT  = (lnQT2 - lnQT0) /(2*dT);
% d2lnqdT = (lnQT2+lnQT0-2*lnQT1)/dT^2;


Cpint = 2*T*dlnqdT + T^2*d2lnqdT;

cp = Cpint + 2.5;

end
