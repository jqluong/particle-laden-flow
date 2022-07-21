A = set_constants2(30);

phi0 = 0.3;
X0 = .25;
p0 = [0.4 -6];

sol = solve_bidensity_ODE(phi0,X0,p0,A);
Z = sol.Z;
X = sol.X;
phi = sol.phi;
sigma = sol.sigma;
T = sol.T;
figure(1)
clf
plot(Z,phi,'-k',Z,X,'-r');
Zt = Z(Z < T);
Xt = X(Z < T);
figure(2)
clf
%subplot(1,2,1)
%plot(Z,log(X),'-k');
%subplot(1,2,2)
plot(T-Zt,(T-Zt).*log(1-Xt),'.-k');
