% Solves the bidensity IVP with log X
% Input:
%   Y0       -  The initial conditions as (phi(0),log(X(0),sigma(0)) 
%   phi0, X0 -  The given total/lighter concentration
%   A        -  Struct of relevant parameters

% Output:
%   G             -  The goal function for shooting
%   Z             -  The values of s
%   X, phi, sigma -  The output values (as X, not log(X))

%TODO: Merge this with fwd_shoot

function [G, Z,X, phi,sigma] = fwd_solve(Y0,phi0,X0,A)
%Solves the IVP (forward)

ffL = @(z,Y) bidensity_FL(z,Y,A);
s = Y0(1); t = Y0(2);
[Z,Y] = A.ODEsolve(ffL,[0 1],[s t 1+(A.rhos+ (A.rhos-A.rhos)*X0)*phi0],A.vopt);
phi = Y(:,1); sigma = Y(:,3);
X = Y(:,2);

T1 = Z(end);
G = [sigma(end) - (1-T1); trapz(Z,X.*phi)/phi0 - X0];

%Heuristic fix for shooting
%if(phi(end) > A.phimax/2 && T1 < 1)
%    G(2) = G(2) + A.phimax*(1-T1);
%end