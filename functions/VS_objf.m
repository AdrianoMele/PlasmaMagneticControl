function J = VS_objf(k,VSsys,x0,Tspan)
%  J = VS_objf(k,VSsys,x0,Tspan)
%    Objective function for VS optimization procedure

ksys = ss(k,'inputname',{'(IVSU-IVSL)/2','Zdot'},'outputname',{'VSU-VSL'});
VS_clsys = feedback(VSsys,ksys,'name');

y = lsim(VS_clsys,Tspan*0,Tspan,x0,'zoh');

J = 0.5*y(:,1)'*y(:,1) + ... % Imbalance Current
    0*  y(:,2)'*y(:,2) + ... % Zp (unused)
    1*  y(:,3)'*y(:,3);      % Zdot
end

