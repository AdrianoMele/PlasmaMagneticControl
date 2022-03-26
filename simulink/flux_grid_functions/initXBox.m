function [RBox, ZBox, irow, icol, mask] = initXBox(xpGuess, R_grid, Z_grid, NPmax)
% [RBox, ZBox, mask] = initXBox(xpGuess, R_grid, Z_grid, [NPmax])
% xpGuess is a first guess of the X-point position; R_grid and Z_grid are defined in initFluxGrid.
% The returned variables are:
% RBox,ZBox: the grid containing the sensors for the interpolation
% icol, irow: indexes of the selected points on the grid
% mask is a matrix whose elements are 1 if the corresponding element of the
%   grid is in RBox, ZBox, NaN otherwise.

% if nargin < 4 
%     NPmin = 3;
% end
if nargin < 4
    NPmax = 3;
end

Nr = size(R_grid, 2);
Nz = size(R_grid, 1);

[~,iR] = min(abs(R_grid(1,:)-xpGuess(1)));
[~,iZ] = min(abs(Z_grid(:,1)-xpGuess(2)));
np = min([iR-1, iZ-1, Nr-iR, Nz-iZ, NPmax]);

irow = iZ-np:iZ+np;
icol = iR-np:iR+np;

mask = NaN(Nz,Nr);
for i = icol
    for j = irow
        mask(j,i) = 1;
    end
end

[RBox, ZBox] = meshgrid(R_grid(1,icol), Z_grid(irow,1));