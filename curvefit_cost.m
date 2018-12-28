% Function to evaluate curve fit cost function
function Error = curvefit_cost(X,T_given)

global data

Error = 0;

for T = T_given
    
    Cp     = spline(data(:,1),data(:,2),T);
    Cp_fit = X(1)/T^2 + X(2)/T + X(3) + X(4)*T + X(5)*T^2 + X(6)*T^3 + X(7)*T^4; 
    Error = Error + (Cp-Cp_fit)^2;
    
end



