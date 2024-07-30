% Computes the bidensity ODE function, 
% except using sLX = log(X)/gamma in place of X
% Here y = (phi,sLX,sigma)
%
% A is the struct containing relevant parameters

function dF= bidensity_FL(~,y,A) %first parameter is time

phimax = A.phimax; gamma = A.gamma;
phi= y(1);
X = y(2);
%X = exp(y(2)*gamma);
sigma= y(3);
d = 3; %Physical Dimension of system

rhoX = A.rhos;

%Auxillary functions
b = abs(A.d1 - A.d2 / (A.d1 + A.d2));
phi_m = A.phimax * (1 + 3/2 * b^(3/2) * (X)^(3/2) * (1 - X));
%mu = (1 - phi/phi_m)^-2;
phi_tr = 0.4;
D_tr = 1/2*min(phi^2,phi_tr^2);
Kc = A.Kc;
Kv = A.Kv;
tol = 1e-4; %tol on how close to 0/1 we're willing to get

if(abs(phi) <= tol || phi <= 0) %Degenerate case,phi >= phimax || 
    dphi = 0;
    dsLX = 0;
elseif (abs(X - 0) < tol || abs(X - 1) < tol) %X = 0,1 solve ODE in 1 species 
    dphi = 1/sigma * (1 + 2*(Kv - Kc)/Kc * phi/(phi_m - phi))^-1 ...
        *( (1 + rhoX*phi)*phi - 2*rhoX*cot(A.alpha)/(9*Kc) * (1-phi));
    dsLX = 0;
else 
    %There's a bunch of terms to build up to create the ODEs
    %Diffusion type coefficient matrices
    %Construct AA, matrix for drift
    %Implicitly divide out a d1^2 in the first row and a d2^2 in the second
    AA = zeros(2,2);
    AA(1,1) = 1;
    AA(2,2) = 1;
    AA(1,2) = 1/4 * (A.d1 + A.d2)^2 / 2^(d+1) * (1 + A.d1/A.d2)^d / ...
        (1 + (A.d1/A.d2)^d) * 1/(A.d1^2); 
    AA(2,1) = 1/4 * (A.d1 + A.d2)^2 / 2^(d+1) * (1 + A.d1/A.d2)^d / ...
        (1 + (A.d1/A.d2)^d) * 1/(A.d2^2); 
    %Construct G, matrix to invert
    xi = phi/phi_m * (phi_m - A.phimax)*(3 - 5*X)/(2*X*(1-X));
    g = Kc*(sigma - phi*sigma*2*(phi_m-phi)^-1) + ...
        Kv*sigma*phi*2*(phi_m - phi)^-1;
    h1 = Kc*(sigma*phi + phi*X*sigma*2*xi*(phi_m - phi)^-1); 
    h2 = Kc*(-sigma*phi + phi*(1-X)*sigma*2*xi*(phi_m - phi)^-1); 
    G = -[X*g h1 ; (1-X)*g h2];
    L = -1/4*sigma*D_tr*phi*[0 1;0 -1]; %Represents tracer
    %Form ODE
    grav = 2*cot(A.alpha)/9 * (1 - phi)*rhoX*[phi*X; (1-X)*phi];
    v = ( diag([phi*X, (1-phi)*X])*AA*G+L)^-1 * (grav ...
        + diag([phi*X, (1-phi)*X])*AA*[phi*X*(-1-phi*rhoX) ; phi*(1-X)*(-1-phi*rhoX)] ); 
    
    
    %Time to actually make the ODEs
    dphi = v(1);
    
    if(X<=0 || X>=1) 
        dsLX = 0;
    else
       %dsLX = 1/gamma*1/X*(L(phi,X,sigma)*Q(phi,X) - P(phi,X)*N(phi,X,sigma)) ...
        %/(-M(phi,X,sigma)*N(phi,X,sigma) + L(phi,X,sigma)*O(phi,X,sigma));
        dsLX = v(2);
    end
end

dsigma= -1-rhoX*phi;

dF=[dphi dsLX dsigma]';
end

%Auxillary functions
