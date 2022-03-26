function [nlc, nlceq] = VS_nlcon(k,VSsys)
% [nlc,nlceq] = VS_nlcon(k,VSsys)
%   Nonlinear constraint on maximum system eigenvalues

ksys = ss(k,'inputname',{'(IVSU-IVSL)/2','Zdot'},'outputname',{'VSU-VSL'});
VS_clsys = feedback(VSsys,ksys,'name');

% continuous-time system eigenvalues
nlc = max(real(eig(VS_clsys.A))); % each superconductive coil will give an eigenvalue in 0
nlceq = [];

if isempty(nlc)
  keyboard
end
end

