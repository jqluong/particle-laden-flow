function phimax = calc_phimax(phi1,phi2,A)
    d1 = A.d1;
    d2 = A.d2;
    b = (d1 - d2)/(d1 + d2);
    phimax0 = A.phimax;
    phi = phi1 + phi2;
    phimax = phimax0*(1 + 3/2 * b^(3/2)*(phi1./phi).^(3/2).*(phi2./phi));
end