% Builds a flux table for interpolation, with x-coordinate phi0
% and y-coordinate X0, and fluxes F(:,:,1:3). 
%
% Varargin stores densities, etc. passed to set_constants (if something
% other than the default is needed). 
%
% N is the number of grid points in each direction.
%
% WARNING: Might fail for very large phi. May be inaccurate there.
% Also note: IF the density difference is large, the solver becomes sad.
% This has not yet been fixed.

function [phi0, X0, F, failed_pts] = build_fluxtableF(ang,N,varargin)

A = set_constants2(ang,varargin{:});

%Adjustable solver settings
A.vopt.RelTol = 1e-7;
A.vopt.AbsTol = 1e-7;
A.ODEsolve = @ode15s;
A.max_it = 60;

phi0 = linspace(0,A.phimax,N);
X0 = linspace(0,1,N);

F = zeros(N,N,3);
g = [0.5 0.5];
failed_pts = [];

%Trivial cases
F(:,1,1) = 1/3;
F(:,1,2:3) = 0;
F(:,end,:) = 0;

%gamma is a scaling factor to make shooting go more smoothly; we set 
%Y = log(X)/gamma and solve the system with Y. IT is chosen so that
%Y(0) is approximately 1 for a solution to make the shooting solver happy. 
gamma0 = A.gamma; 

%iterate from small phi -> large phi
for j=2:N-1
    fprintf('\n****At %f *****\n',phi0(j));
    
    if(j==2)
        A.gamma = gamma0; % generic gamma
    else
        A.gamma = 1.3*gam_old; %Use the gamma from previous phi, but adjusted
                               % a bit upwards
    end
    
    %Start with the X0 = 1 case
    fprintf('***At %f %f\n',phi0(j),X0(end));
    g = [g(1) 0];
    sol = solve_bidensity_ODE(phi0(j),X0(end),g,A);
    F(N,j,:) = compute_fluxes2(sol);
    
    g = sol.g; 
    
    g(2) = -0.05; %Modify initial guess now that X(0) isn't zero
    
    for i=N-1:-1:1 %Work backwards from X0 = 1 to X0 = 0
        
        if(i < N-1)
           A.gamma = -g(2)*A.gamma; %update so that g(2) is approx. -1
        end
        
        if(i==N-2)
            gam_old = 1.1*A.gamma;
        end
        
        fprintf('***At %f %f\n',phi0(j),X0(i));
        
        try
            sol = solve_bidensity_ODE(phi0(j),X0(i),g,A);
        catch
            fprintf('Failed to solve.\n');
            sol.found = 0;
            %sol.g = g;
        end
        
        if(~sol.found)
            failed_pts = [failed_pts; phi0(j), X0(i)]; %#ok<AGROW>
            g = [0.5 -5];
        else
            g = sol.g; %update guess
        end
        
        fprintf('g: %f\t%f, it = %d\n',g(1),g(2),sol.it);
        
        F(i,j,:) = compute_fluxes2(sol);
    end
    
    %save data into a temp file in case it is needed
    save('temp_data.mat','failed_pts','phi0','X0','F','i');
end