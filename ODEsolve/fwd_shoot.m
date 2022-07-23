%see fwd_solve
%
%TODO: merge with fwd_solve

function G = fwd_shoot(Y0,phi0,X0,A)

ffL = @(z,Y) bidensity_FL(z,Y,A);
s = Y0(1); t = Y0(2);

[Z,Y] = A.ODEsolve(ffL,[0 1],[s t 1+(A.rhos+ (A.rhos-A.rhos)*X0)*phi0],A.vopt);

phi= Y(:,1); 
%X = exp(Y(:,2)*A.gamma);
sigma = Y(:,3);

T1 = Z(end); 
G = [sigma(end) - (1-T1); trapz(Z,phi) - phi0];

% A heuristic cheat if the ODE becomes singular before reaching s=1;
% Assume that phi = phimax and estimate the value of sigma(1)

if(phi(end) > A.phimax/2 && T1 < 1)
    G(2) = G(2) + A.phimax*(1-T1);
end