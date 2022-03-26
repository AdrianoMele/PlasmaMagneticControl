function [v, r_grid, z_grid] = grid2vector(V, R_grid, Z_grid)
% V = grid2vector(V, r_grid, z_grid)
% Takes the matrix V and transforms it into a vector that takes values over
% r_grid and z_grid

Nr = size(R_grid,2);
Nz = size(Z_grid,1);
if(size(V,1)*size(V,2)) ~= Nr*Nz, error('Wrong size'); end

v = zeros(Nr*Nz,1);
r_grid = zeros(Nr*Nz,1);
z_grid = zeros(Nr*Nz,1);
for i = 1 : Nr
    for j = 1 : Nz
        v((i-1)*Nr+j)= V(j,i);
        r_grid((i-1)*Nr+j)= R_grid(j,i);
        z_grid((i-1)*Nr+j)= Z_grid(j,i);
    end
end