% Initialize grid
[r_grid, z_grid, R_grid, Z_grid, F0, idxGrid] = initFluxGrid(Input_struct,y_np,y_type, Input_struct.names_sensors, 0);


% Plot machine structure
close all
h = pdemesh(Input_struct.p, Input_struct.e, []);
set(h, 'Color', [.7 .7 .7]);
axis equal
hold on


% Refine grid
r_finer = linspace(r_grid(1),r_grid(end), 100);
z_finer = linspace(z_grid(1),z_grid(end), 100);
[R_finer, Z_finer] = meshgrid(r_finer, z_finer);
FF = interp2(R_grid, Z_grid, F0, R_finer, Z_finer);


% Find separatrix
try
    psib = LinearModel.YEquil(strncmp(LinearModel.OutputsInfo.Name, 'CV-FLUX', 7));
catch
    psib = y_np(strncmp(y_type(:,1), 'psb_c', 5));
end

[c,~] = contour(R_finer, Z_finer, FF, psib*2*pi*[1 1]);
ii = find(c(1,:)<min(r_grid) | c(1,:)>max(r_grid) | c(2,:)<min(z_grid) | c(2,:)>max(z_grid));
c = c(:, setdiff(1:size(c,2),ii));
plot(c(1,:),c(2,:),'sr')

% Compute magnetic field
dr = diff(R_finer(1,:));
dz = diff(Z_finer(:,1));
[FR,FZ] = gradient(FF,dr(1),dz(1));
% From Grad-Shafranov eqn., we know that 
% Bz = (1/2*pi*r)*(dpsi/dr)   and    Br = -(1/2*pi*r)*(dpsi/dz)
Bz =  (1/2/pi)*(FR./R_finer);
Br = -(1/2/pi)*(FZ./R_finer);


% Find null points (imposing grad(psi) = 0)
tol = min([max(max(abs(Br))), max(max(abs(Bz)))])*0.1; % 3 percent of the maximum value
ix = find(abs(Bz)<1e-2 & abs(Br)< tol);
while (length(ix)>2)
    tol = tol*0.9;
    ix = find(abs(Bz)<1e-2 & abs(Br)< tol);
end
rx = R_finer(ix);
zx = Z_finer(ix);
plot(rx,zx, 'xk','MarkerSize', 10);


% Find magnetic axis
imax = find(FF == max(max(FF)));
ra = R_finer(imax);
za = Z_finer(imax);
plot(ra, za, 'xk', 'MarkerSize', 10);


% Find angular coefficients of the cutting lines and select points
for i = 1:length(rx)
    m(i) = -1/((za-zx(i))/(ra-rx(i)));
    plot(R_grid(1,:), m(i)*(R_grid(1,:) - rx(i)) + zx(i), 'k');
    
    if zx(i) > 0
        ii = (c(2,:) <= m(i)*(c(1,:) - rx(i)) + zx(i));
    elseif zx(i) < 0
        ii = (c(2,:) >= m(i)*(c(1,:) - rx(i)) + zx(i));
    end
    c = c(:, ii);
end
% plot(c(1,:),c(2,:),'sc')
patch(c(1,:),c(2,:),'b');



% Find area and volume of the plasma
in = inpolygon(R_finer,Z_finer,c(1,:),c(2,:));
% plot(R_finer(in),Z_finer(in),'sk')
r0 = mean(R_finer(in));
z0 = mean(Z_finer(in));
plot(r0,z0,'sk')
Apl = polyarea(c(1,:),c(2,:));
Vpl = Apl*2*pi*r0;




