%Computes the single-species ODE function (X0=0 or X0=1)
%
%
% phi_shift is used for the integration of the second transient phase
% in the bidisperse PDE simulation (it represents the 'concentration' of
% particles in what is assumed to be the fluid).
%
% Note: this is subject to change (just leave it as zero if not used).

function dF= bidensity_F1(~,y,X0,A)

phimax = A.phimax; 
phi= y(1);
sigma= y(2);

rhoX = A.rhos;%A.rhos2*(1-X0) + A.rhos1*X0;

if(phi >= phimax || phi<= 0)
    dphi = 0;
else 
    dphi=phi+rhoX*(phi.^2 - A.c2*(1-phi-A.phi_shift));
    dphi = dphi/(sigma.*(1+A.c1*(phi./(phimax-(phi+A.phi_shift)))) );
end

dsigma= -1-rhoX*phi;

dF=[dphi dsigma]';