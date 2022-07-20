% Builds a table of fluxes for the 1-d case using
% the framework of the 2-species code (i.e. with forward shooting).

function [phi0, F, failed_pts] = build_fluxtable1(A,X,N)

phi0 = linspace(0,A.phimax,N);

F = zeros(N,2);
g = [0.1 X]; %1d guess
failed_pts = [];

%Trivial cases
F(1,1) = 1/3;
F(1,2) = 0;
F(end,1) = 0;
F(end,2) = 0;

%iterate from small phi -> large phi
for j=2:N-1
    if(mod(j,10)==1)
        fprintf('At phi0 = %f\n',phi0(j));
    end
    
    A.gamma = 1; %ignore this
    
    %This actually does the 1d solve since X = 0 or X = 1
    sol = solve_bidensity_ODE(phi0(j),X,g,A);
    ftmp = compute_fluxes2(sol);
    F(j,1) = ftmp(1);
    if(X==1)
        F(j,2) = ftmp(2);
    else
        F(j,2) = ftmp(3);
    end


    if(~sol.found) %shouldn't happen
        failed_pts = [failed_pts; phi0(j)]; %#ok<AGROW>
        g = [0.5 X];
    else
        g = sol.g; %update guess
    end
        
    
    %save data into a temp file in case it is needed
    %save('temp_data.mat','failed_pts','phi0','F','j');
end