A = set_constants2(45);

phi0 = 0.75;
X0 = 0.65;
p0 = [phi0 X0];  %Initial guess good enough to be phi0 and X0 for now

%Solve ODE
sol = solve_bidensity_ODE(phi0,X0,p0,A);
Z = sol.Z;
X = sol.X;
phi = sol.phi;
sigma = sol.sigma;
T = sol.T;
figure(1)
clf
% plot(Z,phi,'-k',Z,X,'-r');
% title('phi v. Z')
hold on
plot(Z,X .* phi, 'LineWidth', 1.5, 'DisplayName', 'Species 1')
plot(Z, (1-X) .* phi, 'LineWidth', 1.5, 'DisplayName', 'Species 2')
legend
title("\phi_1, \phi_2 v. z for phi_0: " + phi0 + " and \chi_0: " + X0 + " and \alpha: " + A.ang)
hold off
Zt = Z(Z < T);
Xt = X(Z < T);
figure(2)
clf
%subplot(1,2,1)
%plot(Z,log(X),'-k');
%subplot(1,2,2)
plot(T-Zt,(T-Zt).*log(1-Xt),'.-k');
figure(3)
plot(Z,sigma, 'LineWidth', 1.5)
title("\sigma v. z for \phi_0: " + phi0 + " and \chi_0: " + X0 + " and \alpha: " + A.ang)