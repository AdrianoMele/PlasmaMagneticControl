function [r_grid, z_grid, R_grid, Z_grid, F0, idxGrid] = initFluxGrid(Input_struct,y_np,y_type, mag, plotFlag)

if nargin < 5
    plotFlag = 0;
end

% Flux Sensors Grid
% [n, r, z] = EASTmodel.getSensorsInfo;

n = Input_struct.names_sensors;
r = Input_struct.r_sens;
z = Input_struct.z_sens;

r_grid = r(strncmp(n,'Flux_grid_', 10));
z_grid = z(strncmp(n,'Flux_grid_', 10));
clear n r z

equilFluxes = y_np(strncmp(y_type(:,1),'Flux_grid_', 10));

Nr = sqrt(length(r_grid)); % number of points along the z direction
Nz = sqrt(length(z_grid)); % number of points along the r direction
% R_grid = zeros(Nz,Nr);
% Z_grid = R_grid;
% F0 = R_grid;
% for i = 1 : Nr
%     for j = 1 : Nz
%         R_grid(j,i) = r_grid((i-1)*Nr+j); % it should be equal to meshgrid(unique(r_grid), unique(z_grid));
%         Z_grid(j,i) = z_grid((i-1)*Nr+j);
%         F0(j,i) = equilFluxes((i-1)*Nr+j);
%     end
% end
[F0, R_grid, Z_grid] = vector2grid(equilFluxes, r_grid, z_grid);

idxGrid = find(strncmp(mag,'Flux_grid_', 10));


if plotFlag
    figMap = figure('Position', [450 50 450 600]);
    h = pdemesh(Input_struct.p, Input_struct.e, []);
    set(h, 'Color', [.7 .7 .7]);
    axis equal
    hold on
%     plot3(r_grid,z_grid,expGridSensors,'ob');
    plot3(r_grid,z_grid,equilFluxes,'or')
    legend('Mesh','equil')
    % mesh(R_grid, Z_grid, F0);
    % mesh(R_grid, Z_grid, expMap)
end









% Spare lines

% psib = interp1(controlData.Efit.FluxMap.Time, controlData.Efit.FluxMap.PsiB, T);
% psib = psib*2*pi;
% mesh(controlData.Efit.FluxMap.R, controlData.Efit.FluxMap.Z, fluxMap')

% Flux error plot
% mesh(R_grid, Z_grid, expMap-F0);