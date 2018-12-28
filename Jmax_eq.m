%% Function to evaluate right hand side of Jmax.
% This is done by equating the max rotovibrational energy to the max
% vibrational energy.

function [ret] = Jmax_eq(x, D0, b, re, Mu, Evib, Bv, Dv)

global h c

rm=x(1);
J=x(2);

UJ = D0*(1-exp(-b*(rm-re)))^2 + (0.01)*(h*J*(J+1))/(8*pi^2*c*Mu*(rm/100)^2);

ret(1) = 2*D0*b*((exp(-b*(rm-re))) - exp(-2*b*(rm-re))) - (1e-4)*(2*h*J*(J+1))/((rm/100)^3 * 8*pi^2 *c*Mu);

ret(2) = -UJ + Evib + Bv*J*(J+1)- Dv*(J^2)*(J+1)^2;

end