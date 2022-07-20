function F = compute_fluxes2(sol)
% takes a solution struct from the ODE solver solve_bidensity_ODE, 
% and returns a *row* vector of fluxes for the normalized ODE system.

F = zeros(1,3);

T = sol.T;
if(isempty(T))
    T = 1;
end
im = find(sol.Z>=T,1);
if(isempty(im))
    im = length(sol.Z);
end
Z = sol.Z(1:im);
u = sol.u(1:im);
phi = sol.phi(1:im);
X = sol.X(1:im);

F(1) = trapz(Z, u);
F(2) = trapz(Z, X.*phi.*u);
F(3) = trapz(Z, phi.*u) - F(2);

if(T < 1)
    %Do the rest of the first flux analytically
    F(1) = F(1) + u(end)*(1-T) + 1/3*(1-T)^3;
end