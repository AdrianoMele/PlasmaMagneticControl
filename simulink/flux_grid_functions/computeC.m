function [ C_fs, C_br, C_bz ] = computeC( flux_sens, b_sens, Cf, f0, r_grid, z_grid)
% [ C ] = computeC( flux_sens, b_sens, C_flux_sens, f0, r_grid, z_grid)
% Computes the C matrix in a given set of flux and field points
% . flux sens contains the coordinates of the flux sensors (1st column -> r,
%       2nd column -> z)
% . b_sens contains the coordinates of the field sensors
% . Cf is the C matrix associated to the grid of flux sensors
% . f0 is the vector of the equilibrium grid of flux sensors
% . r_grid and z_grid are the coordinates related to the sensors in the grid

if size(flux_sens, 2) ~= 2
    error('Wrong size')
end

% addpath(['..' filesep 'functions'])
[F0, R_grid, Z_grid] = vector2grid(f0, r_grid, z_grid);

nStates = size(Cf, 2);

C_fs = zeros(size(flux_sens,1), nStates);
if isempty(b_sens)
    for i = 1 : nStates
%         dx = zeros(nStates, 1);
%         dx(i) = 1;
        df = Cf(:,i);        
        [dF, ~, ~] = vector2grid(df, r_grid, z_grid);       
        C_fs(:,i) = interp2(R_grid, Z_grid, dF, flux_sens(:,1), flux_sens(:,2), 'linear'); 
        C_br = [];
        C_bz = [];
    end

else
    dr = diff(R_grid(1,1:2));
    dz = diff(Z_grid(1:2,1));
    [FR,FZ] = gradient(F0,dr,dz);
    Bz0 =  (1/2/pi)*(FR./R_grid);
    Br0 = -(1/2/pi)*(FZ./R_grid);
    
    for i = 1 : nStates
        dx = zeros(nStates, 1);
        dx(i) = 1;
        df = Cf*dx;
        
        [FF, ~, ~] = vector2grid(f0+df, r_grid, z_grid);
               
        [FR,FZ] = gradient(FF,dr,dz);
        Bz =  (1/2/pi)*(FR./R_grid);
        Br = -(1/2/pi)*(FZ./R_grid);
        
        C_fs(:,i) = interp2(R_grid, Z_grid, FF-F0, flux_sens(:,1), flux_sens(:,2), 'linear'); 
        C_br(:,i) = interp2(R_grid, Z_grid, Br-Br0, b_sens(:,1),   b_sens(:,2), 'linear'); 
        C_bz(:,i) = interp2(R_grid, Z_grid, Bz-Bz0, b_sens(:,1),   b_sens(:,2), 'linear'); 
               
    end
    
end


% rmpath(['..' filesep 'functions'])
end




