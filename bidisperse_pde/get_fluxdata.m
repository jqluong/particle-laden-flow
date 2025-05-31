% Build the flux data from a table of pre-computed data
% which must be named ftable_[ang].mat
%
% The flux function bF is called with F(Y,fpars{:}) and the 
% Jacobian bJac is called with J(Y,fpars{:},dfpars{:}).
%
% fpars is a cell array containing the grid points phi0_t and X0_t
% and the fluxes as an N by M by 3 array.
% dfpars contains three N by M by 2 arrays; the last dimension specifies
% deriv. with respect to phi or X
% (so dfpars{1} is df with df(:,:,1) = df/dp and df(:,:,2) = df/dx

function [fpars, dfpars] = get_fluxdata(vname,varargin)
%varargin allows a moving average to be specified for the derivative,
%which may (needs to be tested) improve smoothness.
%The code is taken from FEX (#16997, John D'Errico)
%
% Note: Ignore the above and varargin

%load the pre-computed data 
ftable = load(vname + ".mat",vname); ftable = ftable.(vname);
phi0_t = ftable.phi0;
X0_t = ftable.X0;
flux_values = ftable.F; 

fpars = {phi0_t,X0_t,flux_values}; %Parameter list for bF
dphi = phi0_t(2) - phi0_t(1);
dX = X0_t(2) - X0_t(1);

if(nargin > 1) %fix this later
    ord = 3;
    window_length = 7;
    dfp = zeros(length(X0_t),length(phi0_t));
    dg1p = zeros(length(X0_t),length(phi0_t));
    dg2p = zeros(length(X0_t),length(phi0_t));
    dfx = zeros(length(X0_t),length(phi0_t));
    dg1x = zeros(length(X0_t),length(phi0_t));
    dg2x = zeros(length(X0_t),length(phi0_t));
    
    for i=1:length(X0_t)
       dfp(i,:) = movingslope(flux_values(i,:,1),window_length,ord,dphi);
       dg1p(i,:) = movingslope(flux_values(i,:,2),window_length,ord,dphi);
       dg2p(i,:) = movingslope(flux_values(i,:,3),window_length,ord,dphi);
    end
    
    for j=1:length(phi0_t)
       dfx(:,j) = movingslope(flux_values(:,j,1),window_length,ord,dX);
       dg1x(:,j) = movingslope(flux_values(:,j,2),window_length,ord,dX);
       dg2x(:,j) = movingslope(flux_values(:,j,3),window_length,ord,dX);
    end
    
else    
    [dfp, dfx] = gradient(flux_values(:,:,1),dphi,dX);
    [dg1p, dg1x] = gradient(flux_values(:,:,2),dphi,dX);
    [dg2p, dg2x] = gradient(flux_values(:,:,3),dphi,dX);
end

df(:,:,1) = dfp;
df(:,:,2) = dfx;
dg1(:,:,1) = dg1p;
dg1(:,:,2) = dg1x;
dg2(:,:,1) = dg2p;
dg2(:,:,2) = dg2x;

dfpars = {df, dg1, dg2};