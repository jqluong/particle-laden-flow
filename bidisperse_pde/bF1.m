%Flux for the bidisperse problem, takes in Y = (h,n1,n2)
%flux_values has flux_values(i,j,:) = flux at phi0(j), X0(i)

function F = bF1(Y, phi0_t,flux_values_1)

F = zeros(2,1);
if(Y(1) < 1e-12)
    phi = 0;
else
    phi = Y(2)/Y(1);
end

F(1) = interp1(phi0_t,flux_values_1(:,1),phi);
F(2) = interp1(phi0_t,flux_values_1(:,2),phi);
F = Y(1)^3*F;

end