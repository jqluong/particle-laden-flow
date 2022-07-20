% Script to generate a flux table and save it in the proper format.
% Can be edited as needed.

varargs = {}; %arguments to be used by set_constants (densities, etc.)

%fname = ['ftable_',num2str(ang)];
fname = 'ftable_15'; %save the flux using this name for the struct *and* filename
ang = 15; %angle
N = 100; %Number of grid points per dimension

A = set_constants2(ang,varargs{:}); %varargs allows for different densities

[phi0, X0, F, failed_pts] = build_fluxtableF(ang,N,varargs{:});

ftable.phi0 = phi0;
ftable.X0 = X0;
ftable.F = F;
ftable.ang = ang;
ftable.rhos1 = A.rhos1;
ftable.rhos2 = A.rhos2;
ftable.A = A;

% cheat to make sure it has the right name
% (this is probably a terrible way to do it).
assignin('base',fname,ftable); 

save([fname,'.mat'],fname); %save into a file with the same name

