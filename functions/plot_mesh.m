function hm = plot_mesh(Input_struct)
%[hm,hp] = plot_mesh(Input_struct)
%   Plot mesh and return plot handle.
hm = pdemesh(Input_struct.p, Input_struct.e, []);
set(hm,'Color',[.7 .7 .7]);

end

