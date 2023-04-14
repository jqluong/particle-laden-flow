function [val,ist,dir] = phi_event_1(z,y,A)
%stopping criteria for the one-species ODE solve
M = 1 + A.rhos*A.phimax;
ep = 1e-8;
val(1) = +((y(1)>ep)&&(y(1)<A.phimax-ep)); %stop if phi = 0 or phim.
val(2) = +((y(2)>ep)||(z>1-M*ep)); %stop if sigma < ep and z < 1 - M*ep.
ist = [1 1]; dir = [0 0];