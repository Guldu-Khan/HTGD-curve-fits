% Function to evaluate rhs of composition

function X = comp_rhs(x,P0_N2,P0_O2,KpN2,KpO2,KpNO,KpNOe)
% function X = comp_rhs(x,P0_N2,P0_O2,KpN2,KpO2,KpNO)

% X = zeros(7,1);
X = zeros(5,1);

% Species conservation

X(1) = (2*x(1) + x(3) + x(4) + x(5) - 2*P0_N2)^2;
X(2) = (2*x(2) + x(3) + x(4) + x(6) - 2*P0_O2)^2;

% X(1) = 2*x(1) + x(3) + x(4) - 2*P0_N2;
% X(2) = 2*x(2) + x(3) + x(5) - 2*P0_O2;

% Charge conservation
X(3) = (x(4) - x(7))^1;

% Equilibrium constants
X(4) = (KpN2*x(1) - x(5)^2)^2;
X(5) = (KpO2*x(2) - x(6)^2)^2;
X(6) = (KpNO*x(3) - x(5)*x(6))^2;
X(7) = (KpNOe*x(3) - x(4)*x(7))^2;

% X(3) = KpN2*x(1) - x(4)^2;
% X(4) = KpO2*x(2) - x(5)^2;
% X(5) = KpNO*x(3) - x(4)*x(5);

end
