% Computes the bidensity ODE function, 
% except using sLX = log(X)/gamma in place of X
% Here y = (phi,sLX,sigma)
%
% A is the struct containing relevant parameters

function dF= bidensity_FL(~,y,A) %first parameter is time

phimax = A.phimax; gamma = A.gamma;
phi= y(1);
%sLX = y(2);
X = y(2);
%X = exp(y(2)*gamma); %We work with normal X in this ODE, and convert back
sigma= y(3);
d = 3; %Physical Dimension of system

rhoX = A.rhos;

%Auxillary functions and variables
b = abs(A.d1 - A.d2 / (A.d1 + A.d2));
phi_m = phimax * (1 + 3/2 * b^(3/2) * (X)^3/2 * (1 - X));
%mu = (1 - phi/phi_m)^-2; %mu gets cancelled out in equations
phi_tr = 0.4;
D_tr = 1/2*min(phi^2,phi_tr^2);
Kc = A.Kc;
Kv = A.Kv;
tol = 1e-7; %tol on how close to 0/1 we're willing to get
alpha = A.alpha; %Already given in radians

dsigma= -1-rhoX*phi;
if(abs(phi) <= tol || phi <= 0) %Degenerate case ,  phi >= phi_m + tol || 
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
    d1 = A.d1 * 10^8;
    d2 = A.d2 * 10^8;
    a = (d1 + d2)^2 /2^(d+1) * (d1+d2)^d/(d1^d + d2^d);
    detA = d1^2 * d2^2 - a^2;
    %Construct G, matrix to invert
    xi = phi/phi_m * (phi_m - A.phimax)*(3 - 5*X)/(2*X*(1-X));    
    %Time to actually make the ODEs 
    %chi'
    den_X1 = detA*Kc*phi*X*(1-X);
    den_X1 = den_X1 + D_tr*(d1^2*d2^2*(1-X)^2+a*(d1^2+d2^2)*(1-X)*X+d1^2*d2^2*X^2);
    num_X = ((2*X-1)*d1^2*d2^2 + a*(d1^2 - (d1^2 + d2^2)*X));
    dsLX = 8*cot(alpha)/9 * rhoX*((1-phi)*(1-X)*X*num_X) /(sigma*den_X1);
    %dsLX = v(2)/(X*gamma);
    
    %phi'
    g = Kc*sigma*(phi_m - phi) + (Kv - Kc)*2*phi*sigma;
    h1 = Kc*(phi*sigma*(phi_m-phi) + phi*X*sigma*xi);
    h2 = Kc*(-phi*sigma*(phi_m-phi) + phi*(1-X)*sigma*xi);
    den_phi = (X^2 + a*(1/d1^2 + 1/d2^2)*X*(1-X) + (1-X)^2);
    num_phi1 = (X + a/d2^2*(1-X))*h1 + ((1-X) + a/d1^2*X)*h2;
    num_phi2 = 8/9*cot(alpha)*rhoX*(1-phi) + ...
        phi*dsigma*(X^2 + a*(1/d1^2 + 1/d2^2)*X*(1-X) + (1-X)^2);
    dphi = -(num_phi1*dsLX + num_phi2*(phi_m - phi))/(g*den_phi);
end



dF=[dphi dsLX dsigma]';
end

%Auxillary functions
