%%stopping criteria for the bidiserpse ODE
function [val,ist,dir] = phi_event(z,y,A)
M = 1 + A.rhos*A.phimax;
ep = 1e-8;
b = abs(A.d1 - A.d2 / (A.d1 + A.d2));
X = y(2);
phi_m = A.phimax * (1 + 3/2 * b^(3/2) * (X)^3/2 * (1 - X));
val(1) = +((y(1)>ep)); %stop if phi = 0 or phim. AND (y(1)<phi_m-ep
val(2) = 1;
val(3) = +((y(3)>ep)||(z>1-M*ep)); %stop if sigma < ep and z < 1 - M*ep.
ist = [1 0 1]; dir = [0 0 0];