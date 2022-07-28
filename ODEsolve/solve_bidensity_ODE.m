%Solves the equilbrium ODE with initial concentration phi0 and fraction X0.
%The initial guess is p0 = [phi(0) X(0)].

%returns a struct with the solution and a guess for later use (sol.g).
%Failure of the solver is indicated by sol.it returning -1.
%
% Inputs:
%   phi0, X0    - The initial concentration/fraction
%          A    - Struct of ODE parameters from set_constants2
%         p0    - The initial guess for (p(0),X(0)) OR (p(0),log(X(0)).
%                 The choice is processed by the sign of p0(2).
%
% Outputs:
%   sol
%       .Z      - The integration points. Note: the integration stops when 
%                 phi = 0, so Z(end) is not always 1.
%       .X, 
%       .phi, 
%       .sigma   - The ODE solution; sigma is the shear stress.
%       .u       - The velocity, from integrating mu^(-1)*sigma
%       .T       - The first z-value with phi = 0 (so T should be Z(end)).
%       .it      - # of iterations used in shooting.

function sol = solve_bidensity_ODE(phi0,X0,p0,A)



tol = 100*A.vopt.RelTol; %can be changed

%---------------
%Process initial guess...

if(size(p0,1)==1)
    p0 = p0';
elseif(isempty(p0))     
    p0 = [phi0; X0]; 
end

if(p0(1) < 0 || p0(1) > A.phimax)
    p0(1) = A.phimax/2;
end

if(p0(2) > 0) %Change to SLX = log(X)/gamma
    p0(2) = log(p0(2))/A.gamma;
end

sol_found = 0;
%----------------------

if(X0==0 || X0==1)
    % The one-species case:
    
    phi_e = @(~,y) phi_event_1(0,y,A);
    A.vopt.Events = phi_e; %Stopping criterion
    F = @(s) fwd_shoot_1(s,phi0,X0,A);
    [s,it] = bisect(F,0,A.phimax,tol);
    [~,Z,phi,sigma] = fwd_solve_1(s,phi0,X0,A);
    X(1:length(Z),1) = X0;
    g = [s; X0];
    if(it<=A.max_it)
        sol_found = 1;
    else
        fprintf('Solver failed for some reason. Results may be inaccurate.\n');
        sol_found = 0;
    end
else
    % Two species case
    
    phi_e = @(~,y) phi_event(0,y,A);
    A.vopt.Events = phi_e; %Stopping criterion
    %Ensure X-guess is not zero or one
    if(p0(2)==0 || p0(2)==1 || p0(2) == -Inf)
          p0(2) = -1;
    end
    
    %F = @(Y) fwd_shoot(Y,phi0,X0,A);
    pars = {phi0,X0,A}; %paramters for the shooting function
    bounds = [0 A.phimax; -Inf 0];

    p = p0;
    k = 0; max_k = 20;
    while(~sol_found  && k < max_k)
        k = k + 1;
        [Y,it] = nlsolve_broyden(p,@euler,2,tol,bounds,A.max_it,pars);
        if(it==-1 || it>=A.max_it) %solver failure
            %try a new guess?
            if(k==max_k) %this isn't attempted
                p = [0.5 log(0.5)]';
            else
                xi = 0.05;
                p(1) = p(1) + xi*(A.phimax - p(1));
                p(2) = max(0.9*p(2),-50);
            end
                %...could try something else here
        else 
            sol_found = 1; 
        end
    end
  
    if(k>=max_k)
        fprintf('Warning: solver failed for some reason. Results may be inaccurate.\n');
    end
    
    [~,Z,X,phi,sigma] = fwd_solve(Y,phi0,X0,A);
    g = Y; %Shooting value that solves the problem
end

%--------------
%Build the solution structure
sol.it = it; sol.Z = Z;
sol.X = X; sol.phi = phi; sol.sigma = sigma;
sol.g = g;

%velocity u in the region where phi is non-zero
sol.u = cumtrapz(Z,(1-(phi+A.phi_shift)/A.phimax).^2.*sigma); %shift for transient 2 case (TODO: document)
sol.found = sol_found;

%Find the extent of the particle layer
%This is important for the flux calculation!
sol.T = sol.Z(find(sol.phi<=1.2e-6,1));
if(isempty(sol.T))
    sol.T = sol.Z(end);
end