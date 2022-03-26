function [V, R_grid, Z_grid] = vector2grid(v, r_grid, z_grid)
% [V, R_grid, Z_grid] = vector2grid(v, R_grid, Z_grid)
% Takes the vector v and transforms it into a matrix tht takes values over
% R_grid and Z_grid

% if(length(v)) ~= Nr*Nz, error('Wrong size'); end
if (length(v))~=length(r_grid) || (length(v))~=length(z_grid), error('Wrong size'); end

Nr = sqrt(length(r_grid));
Nz = sqrt(length(z_grid));

V = zeros(Nr,Nz);
R_grid = zeros(Nr,Nz);
Z_grid = zeros(Nr,Nz);
for i = 1 : Nr
    for j = 1 : Nz
        V(j,i) = v((i-1)*Nr+j);
        R_grid(j,i) = r_grid((i-1)*Nr+j);
        Z_grid(j,i) = z_grid((i-1)*Nr+j);
    end
end