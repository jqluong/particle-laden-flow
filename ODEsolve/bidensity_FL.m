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

rhoX = A.rhos;

if(phi >= phimax || phi<= 0) %Degenerate case
    dphi = 0;
    dsLX = 0; 
else 
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