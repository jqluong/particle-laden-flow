%This is a script from making a quick gif from a sol file
%Assume you have a sol struct
run_name = "presentation_low_chi";
n = length(sol.t_out);
V_slurry = sol.Vsl*10^6;
g = 9.8;
Vsl = (V_slurry)*(1e-6); %slurry volume in m^3;
Ltr = 1; %track length in m
angr = ang*pi/180; %angle in radians
wtr = (14)*1e-2; % track width in m
H = Vsl/(wtr*Ltr); %char. height H
t_dim = A.nul*Ltr/(H^2*g*sin(angr)); %time scaling (in s)
for i = 1:n
    plot(sol.X,sol.Y(:,:,i),'LineWidth',4);
    xlabel('Track length (m)')
    xlim([-0.1 1.3])
    ylabel('Height profile')
    legend("Fluid", "Species 1 (Larger)", "Species 2 (Smaller)")
    fontsize(40, "points")
    colororder(["#000000" "#20ADFF" "#FF0000"]) %fluid, larger, smaller
    %colororder(["#000000" "#FF0000" "#0000FF"]) %black, red, blue
    title(sprintf('t = %5.2f s',sol.t_out(i)*t_dim));
    set(gca,'TickLength',[0.025 0.025],'LineWidth',2);
    fig = gcf;
    fig.WindowState = "maximized";
    exportgraphics(gca,run_name + ".gif","Append",true)
end