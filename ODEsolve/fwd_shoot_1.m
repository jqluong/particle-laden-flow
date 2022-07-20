% Forward shooting, for phi and sigma in the monodisperse case (X0=0 or X0=1)
% where the densities are specified in the struct A.
%
%goal: satisfy shear-stress constraint sigma(1) = 0

function G = fwd_shoot_1(s,phi0,X0,A)

ff1 = @(z,Y) bidensity_F1(z,Y,X0,A);

[Z,Y] = A.ODEsolve(ff1,[0 1],[s 1+(A.rhos2+ (A.rhos1-A.rhos2)*X0)*phi0],A.vopt);

sigma = Y(:,2);

T1 = Z(end); 
G = sigma(end) - (1-T1);