function hp = plot_plasma(Input_struct,psi,psicont)
%hp = plot_plasma(Input_struct,psi,psicont)
%   Plot plasma contours and return plot handle(s). 
hold on
hp = pdecont(Input_struct.p, Input_struct.t, psi, psicont);
end

