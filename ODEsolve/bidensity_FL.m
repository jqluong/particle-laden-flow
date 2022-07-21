% Computes the bidensity ODE function, 
% except using sLX = log(X)/gamma in place of X
% Here y = (phi,sLX,sigma)
%
% A is the struct containing relevant parameters

function dF= bidensity_FL(~,y,A) %first parameter is time

phimax = A.phimax; gamma = A.gamma;
phi= y(1);
%sLX = y(2);
X = exp(y(2)*gamma);
sigma= y(3);
d = 3; %Physical Dimension of system

rhoX = A.rhos;

%Auxillary functions
b = abs(A.d1 - A.d2 / (A.d1 + A.d2));
delta = @(phi,x) A.rho * phi * (x * A.d1/2 + (1-x)*A.d2/2)^2 / ...
    ( (x * A.d1/2)^2 + ( ( 1-x) * A.d2/2)^2);
phi_m = @(phi,x) A.phimax * (1 + 3/2 * b^(3/2) * (x)^3/2 * (1 - x));
S = @(phi,x) (phi_m(phi,x)^2 - phi^2) / (phi_m(phi,x) - phi*(1 - delta));

if(phi >= phimax || phi<= 0) %Degenerate case
    dphi = 0;
    dsLX = 0; 
else 
    %There's a bunch of terms to build up to create the ODEs
    %Diffusion type coefficient matrices
    Kd = A.Kc * 0.61;
    %Construct AA, matrix for drift
    AA(2,2) = zeros(2,2);
    AA(1,1) = 1/4*A.d1^2;
    AA(2,2) = 1/4*A.d2^2;
    AA(1,2) = 1/4 * (A.d1 + A.d2)^2 / 2^(d+1) * (1 + A.di/A.dj)^d / ...
        (1 + (A.di/A.dj)^d);
    AA(2,1) = AA(1,2);
    %Construct DD, matrix for tracer.  These are functions in phi and x
    kb = 1.380649 * (10^-23); %Boltzman constant in m^2 kg / s^2 K
    D1 = kb * 293 / (3 * pi * A.nul * A.rholiq*A.d1);%Assume room temp = 293K
    D2 = kb * 293 / (3 * pi * A.nul * A.rholiq*A.d2);%Assume room temp = 293K 
    DD11 = @(phi,x) Kd * D1 * 1/S(phi,x);
    DD22 = @(phi,x) Kd * D2 * 1/S(phi,x);
    DD12 = @(phi,x) DD11 * phi * x / (phi_m(phi,x) - phi*(1-x))^2;
    DD21 = @(phi,x) DD22 * phi * (1-x) / (phi_m(phi,x) - phi*x)^2;
    
    dphi=phi+rhoX*(phi.^2 - A.c2*(1-phi));
    dphi = dphi/(sigma.*(1+A.c1*(phi./(phimax-phi))) );
    
    if(X<=0 || X>=1) 
        dsLX = 0;
    else
        dsLX = (A.c0/gamma)*(1-X)*phimax/((phimax-phi)*sigma*A.Dtr(phi));
    end
end

dsigma= -1-rhoX*phi;

dF=[dphi dsLX dsigma]';
end

%Auxillary functions
