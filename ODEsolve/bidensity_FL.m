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
phi_m = A.phimax * (1 + 3/2 * b^(3/2) * (X)^3/2 * (1 - X));
mu = (1 - phi/phi_m)^-2; %mu gets cancelled out in equations
phi_tr = 0.4;
D_tr = 1/2*min(phi^2,phi_tr^2);
Kc = A.Kc;
Kv = A.Kv;
tol = 1e-4; %tol on how close to 0/1 we're willing to get
alpha = A.alpha;

if(phi >= phi_m || abs(phi) <= tol || phi <= 0) %Degenerate case,
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
    %Construct G, matrix to invert
    xi = phi/phi_m * (phi_m - A.phimax*(3 - 5*X))/(2*X*(1-X));    
    %Time to actually make the ODEs
    num_1 = (-16/9*cot(alpha)*rhoX*(1 - phi) + phi*(1+phi*rhoX))*(phi_m - phi);
    num_2 = 2*xi*(8/9*cot(alpha)*rhoX*(1-phi)*(2*X-1))/(Kc+(D_tr/phi));
    den = Kc*sigma*(phi_m - phi) + 2*phi*sigma*(Kv - Kv);
    dphi = (num_1 + num_2)/(den);
    
    if(X<=0 || X>=1) 
        dsLX = 0;
    else
        num_X = 8/9*cot(alpha)*rhoX*(1-phi)*(2*X - 1);
        den_X = sigma*(Kc*phi + D_tr);
        dsLX = num_X/den_X;
        %dsLX = v(2)/(X*gamma);
    end
end

dsigma= -1-rhoX*phi;

dF=[dphi dsLX dsigma]';
end

%Auxillary functions
