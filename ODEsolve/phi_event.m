function [val,ist,dir] = phi_event(z,y,A)
%stopping criteria for the bidensity ODE
M = 1 + A.rhos2*A.phimax;
ep = 1e-8;
val(1) = +((y(1)>ep)&&(y(1)<A.phimax-ep)); %stop if phi = 0 or phim.
val(2) = 1;
val(3) = +((y(3)>ep)||(z>1-M*ep)); %stop if sigma < ep and z < 1 - M*ep.
ist = [1 0 1]; dir = [0 0 0];