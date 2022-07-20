%Flux for the bidisperse problem, takes in Y = (h,n1,n2)
%flux_values has flux_values(i,j,:) = flux at phi0(j), X0(i)

function F = bF(Y, phi0_t,X0_t,flux_values)

F = zeros(3,1);
%v = get_p(Y);

v(1) = (Y(2) + Y(3))/Y(1);
v(2) = Y(2)/(Y(2)+Y(3));

F(1) = interp2(phi0_t,X0_t,flux_values(:,:,1),v(1),v(2));
F(2) = interp2(phi0_t,X0_t,flux_values(:,:,2),v(1),v(2));
F(3) = interp2(phi0_t,X0_t,flux_values(:,:,3),v(1),v(2));

F = Y(1)^3*F;

end