function pde_plotter(X,Y,t,t_scale,t_scale_name)
%PDE_PLOTTER Summary of this function goes here
%   Detailed explanation goes here
    plot(X,Y(:,1),'-k',X,Y(:,2),'-r',X,Y(:,3),'-b');
    legend('h','n1','n2','Location','Northwest');
    title(sprintf('t = %f %s',t*t_scale,t_scale_name));
    drawnow
end

