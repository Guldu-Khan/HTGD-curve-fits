%% Function to find the partition function for diatomic gases
function Qint = Q2(data,T)

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

% This is a triple summation of 1.electronic, 2.vibrational and 3.
% rotational partition functions. The rotational and vibrational states are
% coupled and thus there are only two summations - 1.electronic and
% 2.rotovibrational.

% N_max = min(ceil(0.4*T^0.25),size(Te,1));   % For O2
Qint = 0;
% N_max = 15;            % For CO
% N_max = 5;
N_max = size(Te,1);     % N2
% N_max = min(ceil(0.35*T^0.25),size(Te,1));   % For NO
% disp(N_max);

% First summation - electronic
for n = 1:N_max
    
    Qel = g(n)*exp(-(Te(n))*h*c*100/(k*T));
    
    % Cut-off for rotational states - V_max
    E0 = 0.5*We(n) - 0.25*WeXe(n) + 0.125*WeYe(n) + WeZe(n)/16;
    W0 = We(n) - WeXe(n) + 0.75*WeYe(n) + 0.125*WeZe(n);
    W0X0 = WeXe(n) - 1.5*WeYe(n) - 1.5*WeZe(n);
    W0Y0 = WeYe(n) + 2*WeZe(n);
    W0Z0 = WeZe(n);
    b    = We(n)*sqrt(2*pi^2*c*Mu/(100*De(n)*h));

    p = [W0Z0 W0Y0 -W0X0 W0 -(D0(n))];
    r = roots(p);
    r = r(imag(r)==0);
    r = r(r>0);
    V_max = ceil(min(r));
%     disp(V_max)
    
    Qrv = 0;    
    Jguess = 50;                       % Guess value for J_max
%     V_max = 0;
    % Second summation - rotovibrational
    for s = 0:V_max
        
        Evib = E0 + W0*s - W0X0*s^2 + W0Y0*s^3 + W0Z0*s^4;
        Bv   = Be(n) - ae(n)*(s+0.5);
        Dv   = De(n) - beta(n)*(s+0.5);
%         Qvib = exp(-100*(Evib)*(h*c)/(k*T))
        
        % Cut-off for rotational states
        x0      = [1e-3 Jguess];
        options = optimoptions('fsolve','TolFun',1e-9,'MaxFunEvals',1500,'Display','none');
        sol     = fsolve(@Jmax_eq,x0,options,D0(n),b,re(n),Mu,Evib,Bv,Dv);
        J_max   = sol(2);
        
        if J_max < 0
            J_max = 0;
        else
            J_max = ceil(sol(2));
        end
%         disp(J_max)
%         Qrot = 0;
        for J = 0:J_max
            
            Erot = Bv*J*(J+1)-Dv*(J^2)*(J+1)^2;
            Evr  = Evib + Erot;
            Qrv  = Qrv + (2*J+1)*exp(-100*Evr*(h*c)/(k*T));
        end
        
    end
    
    Qint = Qint + Qel*Qrv;
    
end

Qint = Qint/sigma;

end
        
