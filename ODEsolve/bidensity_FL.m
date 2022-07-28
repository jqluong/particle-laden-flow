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
delta = @(phi,x) A.rhos * phi * (x * A.d1/2 + (1-x)*A.d2/2)^2 / ...
    ( (x * A.d1/2)^2 + ( ( 1-x) * A.d2/2)^2);
phi_m = @(phi,x) A.phimax * (1 + 3/2 * b^(3/2) * (x)^3/2 * (1 - x));
S = @(phi,x) (phi_m(phi,x)^2 - phi^2) / (phi_m(phi,x) - phi*(1 - delta(phi,x)));

if(phi >= phimax || phi<= 0) %Degenerate case
    dphi = 0;
    dsLX = 0; 
else 
    %There's a bunch of terms to build up to create the ODEs
    %Diffusion type coefficient matrices
    Kd = A.Kc * 0.61;
    Kc = A.Kc;
    Kv = A.Kv;
    %Construct AA, matrix for drift
    AA = zeros(2,2);
    AA(1,1) = 1/4;
    AA(2,2) = 1/4;
    AA(1,2) = 1/4 * 1/(A.d1^2) * (A.d1 + A.d2)^2 / 2^(d+1) * (1 + A.d1/A.d2)^d / ...
        (1 + (A.d1/A.d2)^d);
    AA(2,1) = 1/4 * 1/(A.d2^2) * (A.d2 + A.d1)^2 / 2^(d+1) * (1 + A.d2/A.d1)^d / ...
        (1 + (A.d2/A.d1)^d);
    %Construct DD, matrix for tracer.  These are functions in phi and x
    kb = 1.380649 * (10^-23); %Boltzman constant in m^2 kg / s^2 K
    D1 = 1/A.d1^2 * kb * 293 / (3 * pi * A.nul * A.rholiq*A.d1);%Assume room temp = 293K
    D2 = 1/A.d2^2 * kb * 293 / (3 * pi * A.nul * A.rholiq*A.d2);%Assume room temp = 293K 
    DD11 = @(phi,x) Kd * D1 * 1/S(phi,x);
    DD22 = @(phi,x) Kd * D2 * 1/S(phi,x);
    DD12 = @(phi,x) DD11 * phi * x / (phi_m(phi,x) - phi*(1-x))^2;
    DD21 = @(phi,x) DD22 * phi * (1-x) / (phi_m(phi,x) - phi*x)^2;
    
    %Build blocks for ODE
    e = @(x,sigma) 2*(Kc - Kv)*sigma*phi;
    xi = @(phi,x,sigma) phi/phi_m(phi,x) * (phi_m(phi,x) - A.phimax*(3 - 5*x))/(2*x*(1-x));
    mu = @(phi,x) (1 - phi/phi_m(phi,x))^2;
    F1 = @(phi,x,sigma) -phi/mu(phi,x) * x * AA(1,1) * (Kc*sigma*x + e(x,sigma)*x) ...
        - phi/mu(phi,x) * x * AA(1,2) * (Kc*sigma*(1-x) + e(x,sigma)*(1-x));
    G1 = @(phi,x,sigma) -phi/mu(phi,x) * x * AA(1,1) * (-Kc*sigma*phi -xi(phi,x,sigma)* e(x,sigma)*x) ...
        - phi/mu(phi,x) * x * AA(1,2) * (-Kc*sigma*phi - xi(phi,x,sigma)*e(x,sigma)*(1-x));
    H1 = @(phi,x) -phi/mu(phi,x) * x * AA(1,1)*(phi*x*(-1 - phi*A.rhos)) ...
        - phi/mu(phi,x)*x*AA(1,2)*(phi*(1-x)*(-1 - phi*A.rhos));
    F2 = @(phi,x,sigma) -phi/mu(phi,x) * (1-x) * AA(1,1) * (Kc*sigma*x + e(x,sigma)*x) ...
        - phi/mu(phi,x) * (1-x) * AA(1,2) * (Kc*sigma*(1-x) + e(x,sigma)*(1-x));
    G2 = @(phi,x,sigma) -phi/mu(phi,x) * (1-x) * AA(1,1) * (-Kc*sigma*phi -xi(phi,x,sigma)* e(x,sigma)*x) ...
        - phi/mu(phi,x) * (1-x) * AA(1,2) * (-Kc*sigma*phi - xi(phi,x,sigma)*e(x,sigma)*(1-x));
    H2 = @(phi,x) -phi/mu(phi,x) * (1-x) * AA(1,1)*(phi*x*(-1 - phi*A.rhos)) ...
        - phi/mu(phi,x)*(1-x)*AA(1,2)*(phi*(1-x)*(-1 - phi*A.rhos));
    J1 = @(phi,x) x*phi / (phi_m(phi,x) - (1-x)*phi^2);
    K1 = @(phi,x) -1 * 2 * cot(A.alpha)*phi*x / 9 * (1 - phi/phi_m(phi,x))*(A.rhos);
    J2 = @(phi,x) (1-x)*phi / (phi_m(phi,x) - x*phi^2);
    K2 = @(phi,x) -1 *2* cot(A.alpha)*phi*(1-x) / 9 * (1 - phi/phi_m(phi,x))*(A.rhos);
    
    %Almost there, these are the big blocks
    L = @(phi,x,sigma) F1(phi,x,sigma) + DD11(phi,x)*x + J1(phi,x)*DD11(phi,x)*(1-x);
    M = @(phi,x,sigma) G1(phi,x,sigma) + DD11(phi,x)*phi + J1(phi,x)*DD11(phi,x)*x;
    P = @(phi,x) -H1(phi,x) - K1(phi,sigma);
    N = @(phi,x,sigma) F2(phi,x,sigma) + DD22(phi,x)*J2(phi,x)*x + DD22(phi,x)*(1-x);
    O = @(phi,x,sigma) G2(phi,x,sigma) + DD22(phi,x)*J2(phi,x)*x - DD22(phi,x)*x;
    Q = @(phi,x) -H2(phi,x) - K2(phi,sigma);
    
    
    %Time to actually make the ODEs
    dphi = (M(phi,X,sigma)*Q(phi,X) - P(phi,X)*O(phi,X,sigma)) ...
        /(M(phi,X,sigma)*N(phi,X,sigma) - L(phi,X,sigma)*O(phi,X,sigma));
    
    if(X<=0 || X>=1) 
        dsLX = 0;
    else
       dsLX = 1/gamma*1/X*(L(phi,X,sigma)*Q(phi,X) - P(phi,X)*N(phi,X,sigma)) ...
        /(-M(phi,X,sigma)*N(phi,X,sigma) + L(phi,X,sigma)*O(phi,X,sigma));
    end
end

dsigma= -1-rhoX*phi;

dF=[dphi dsLX dsigma]';
end

%Auxillary functions
