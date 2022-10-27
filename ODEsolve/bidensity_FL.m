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
phi_m = @(x) A.phimax * (1 + 3/2 * b^(3/2) * (x)^3/2 * (1 - x));
S = @(phi,x) (phi_m(x)^2 - phi^2) / (phi_m(x) - phi*(1 - delta(phi,x)));
mu = @(phi,x) (1 - phi/phi_m(x))^-2;
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
    %Construct G, matrix to invert
    xi = @(phi,x,sigma) phi/phi_m(x) * (phi_m(x) - A.phimax*(3 - 5*x))/(2*x*(1-x));
    g = @(phi,x,sigma) Kc*(sigma - phi*sigma*2*(phi-phi_m(x)^-1)) + ...
        Kv*sigma*phi*2*(phi_m(x) - phi)^-1;
    h1 = @(phi,x,sigma) Kc*(sigma*phi + phi*x*sigma*2*xi(phi,x,sigma)*(phi - phi_m(x))^-1); 
    h2 = @(phi,x,sigma) Kc*(-sigma*phi + phi*(1-x)*sigma*2*xi(phi,x,sigma)*(phi - phi_m(x))^-1); 
    G = @(phi,x,sigma) [x*g(phi,x,sigma) h1(phi,x,sigma) ; (1-x)*g(phi,x,sigma) h2(phi,x,sigma)];
    
    %Form ODE
    grav = @(phi,x) -2*cot(A.alpha)/9 * (1 - phi/phi_m(x))*rhoX*[A.d1^2*phi*x; A.d2^2*(1-phi)*x];
    v = G(phi,X,sigma)^-1 * ( (mu(phi,X)*AA^-1*diag([1/(phi*X), 1/(phi*(1-X))])) * grav(phi,X)...
        + [phi*X*(-1-phi*rhoX) ; phi*(1-X)*(-1-phi*rhoX)] ); 
    %Build blocks for ODE
%     e = @(x,sigma) 2*(Kc - Kv)*sigma*phi;
% 
%     F1 = @(phi,x,sigma) -phi/mu(phi,x) * x * AA(1,1) * (Kc*sigma*x + e(x,sigma)*x) ...
%         - phi/mu(phi,x) * x * AA(1,2) * (Kc*sigma*(1-x) + e(x,sigma)*(1-x));
%     G1 = @(phi,x,sigma) -phi/mu(phi,x) * x * AA(1,1) * (-Kc*sigma*phi -xi(phi,x,sigma)* e(x,sigma)*x) ...
%         - phi/mu(phi,x) * x * AA(1,2) * (-Kc*sigma*phi - xi(phi,x,sigma)*e(x,sigma)*(1-x));
%     H1 = @(phi,x) -phi/mu(phi,x) * x * AA(1,1)*(phi*x*(-1 - phi*A.rhos)) ...
%         - phi/mu(phi,x)*x*AA(1,2)*(phi*(1-x)*(-1 - phi*A.rhos));
%     F2 = @(phi,x,sigma) -phi/mu(phi,x) * (1-x) * AA(1,1) * (Kc*sigma*x + e(x,sigma)*x) ...
%         - phi/mu(phi,x) * (1-x) * AA(1,2) * (Kc*sigma*(1-x) + e(x,sigma)*(1-x));
%     G2 = @(phi,x,sigma) -phi/mu(phi,x) * (1-x) * AA(1,1) * (-Kc*sigma*phi -xi(phi,x,sigma)* e(x,sigma)*x) ...
%         - phi/mu(phi,x) * (1-x) * AA(1,2) * (-Kc*sigma*phi - xi(phi,x,sigma)*e(x,sigma)*(1-x));
%     H2 = @(phi,x) -phi/mu(phi,x) * (1-x) * AA(1,1)*(phi*x*(-1 - phi*A.rhos)) ...
%         - phi/mu(phi,x)*(1-x)*AA(1,2)*(phi*(1-x)*(-1 - phi*A.rhos));
%     J1 = @(phi,x) x*phi / (phi_m(x) - (1-x)*phi)^2;
%     K1 = @(phi,x) -1 * 2 * cot(A.alpha)*phi*x / 9 * (1 - phi/phi_m(x))*(A.rhos);
%     J2 = @(phi,x) (1-x)*phi / (phi_m(x) - x*phi)^2;
%     K2 = @(phi,x) -1 *2* cot(A.alpha)*phi*(1-x) / 9 * (1 - phi/phi_m(x))*(A.rhos);
%     
%     %Almost there, these are the big blocks
%     L = @(phi,x,sigma) F1(phi,x,sigma) + DD11(phi,x,sigma)*x + J1(phi,x)*DD11(phi,x,sigma)*(1-x);
%     M = @(phi,x,sigma) G1(phi,x,sigma) + DD11(phi,x,sigma)*phi + J1(phi,x)*DD11(phi,x,sigma)*x;
%     P = @(phi,x) -H1(phi,x) - K1(phi,sigma);
%     N = @(phi,x,sigma) F2(phi,x,sigma) + DD22(phi,x,sigma)*J2(phi,x)*x + DD22(phi,x,sigma)*(1-x);
%     O = @(phi,x,sigma) G2(phi,x,sigma) + DD22(phi,x,sigma)*J2(phi,x)*x - DD22(phi,x,sigma)*x;
%     Q = @(phi,x) -H2(phi,x) - K2(phi,sigma);
    
    
    %Time to actually make the ODEs
    dphi = v(1);
    
    if(X<=0 || X>=1) 
        dsLX = 0;
    else
       %dsLX = 1/gamma*1/X*(L(phi,X,sigma)*Q(phi,X) - P(phi,X)*N(phi,X,sigma)) ...
        %/(-M(phi,X,sigma)*N(phi,X,sigma) + L(phi,X,sigma)*O(phi,X,sigma));
        dsLX = v(2)/(X*gamma);
    end
end

dsigma= -1-rhoX*phi;

dF=[dphi dsLX dsigma]';
end

%Auxillary functions
