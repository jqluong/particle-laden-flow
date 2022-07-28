%instead of using a matlab default solver im going to hand implement an
%euler's method to debug and see what's happening
%this replaces fwd_shoot
function G = euler(Y0, phi0, X0, A)
    ffL = @(z,Y) bidensity_FL(z,Y,A);
    s = Y0(1); t = Y0(2);
    %Perform an Euler method from 0 to 1
    Z = [0:A.dz:1];
    nodes = length(Z);
    Y = zeros(nodes,3);
    Y(1,1) = s;
    Y(1,2) = t;
    Y(1,3) = 1+A.rhos*phi0;
    for i = 2:nodes
        Y(i,:) = Y(i-1,:) + A.dz*ffL(Z(i-1),Y(i-1,:))';
    end
    %Post processing
    phi= Y(:,1); 
    sigma = Y(:,3);
    T1 = Z(end); 
    G = [sigma(end) - (1-T1); trapz(Z,phi) - phi0];

    % A heuristic cheat if the ODE becomes singular before reaching s=1;
    % Assume that phi = phimax and estimate the value of sigma(1)

    if(phi(end) > A.phimax/2 && T1 < 1)
        G(2) = G(2) + A.phimax*(1-T1);
    end
end