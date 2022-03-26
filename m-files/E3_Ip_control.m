%% Plasma current controller
% In this exercise we want to design a plasma current controller. The $I_p $ 
% control is in charge of keeping the plasma current close to the desired values 
% despite the action of external disturbances (such as additional heating systems) 
% and model uncertainties (e.g. a wrong estimate of the plasma resistivity). Of 
% course, the scenario currents are designed so as to keep $I_p $ close to its 
% reference value during the shot, but in general all of the effects above are 
% only approximately known in advance and hence a closed loop control action is 
% needed.
%% Preliminaries and transformer currents pattern
% As usual, we start by loading our plasma linearized model

addpath ./functions
addpath ./models
addpath ./data 

modelName = fullfile(pwd,'models','SOF@18d66s_FG.mat');
model = load(modelName);

lm.A = -model.L\model.R;
lm.B =  model.L\eye(size(model.L));
lm.C =  model.C;
lm.D =  zeros(size(lm.C,1), size(lm.B,2));
lm.E = -model.L\model.LE;
lm.F =  model.F;
%% 
% Now, the key idea is to find a combination of PF currents capable of acting 
% on the plasma current without affecting the plasma shape too much. 
% 
% *(Notice: for a diverted configuration, the plasma boundary is defined as 
% the isoflux line passing through the X-point. )*
% 
% There are a few different ways to find this combination. The easiest one is 
% to take a few plasma-wall gaps and find a combination of currents which acts 
% on $I_p $ without changing the gaps.
% 
% First we need to quantify the effect of the active currents on the plasma 
% current. We could take the row of the C matrix associated to the plasma current

i_Ip = get_y_idx(model.y_type,'Ipl',1);
C_Ip = lm.C(i_Ip,:);
%% 
% But since Ip is a state variable of our model, this matrix is just a bunch 
% of zeros

C_Ip(:,1:10)
%% 
% So, we need to find an alternative solution. Let's start from the circuit 
% equations
% 
% $$L\dot{I}(t)+RI(t) = V(t)$$
% 
% We are only looking for the effect of the PF currents on $I_p$, so we'll make 
% the following simplifying assumptions:
%% 
% * neglect the effect of the resistivity on both the plasma (very hot) and 
% the coils (superconductive), so $R_{PF} = 0, R_{Ip} \approx 0$
% * neglect the effect of the eddy currents
%% 
% Our circuit equations become
% 
% $$L_{PF}\dot{I}_{PF}(t)+M_{I_p2PF}\dot{I}_p(t)=V_{PF}(t) \\M_{PF2I_p} \dot{I}_p(t)(t)+L_{I_p}\dot{I}_p(t)=0$$
% 
% where the pedices are used to distinguish between the PF currents and the 
% plasma current and the self and mutual inductance term are explicitly separated. 
% 
% From the second equation
% 
% $$dI_p \approx -\frac{M_{PF2I_p}}{L_{I_p} }dI_{PF}(t)$$
% 
% This is (approximately) the contribution of a variation of the PF currents 
% on the plasma current we were looking for! So the matrix we need is

C_PF2Ip   = -model.L(end,1:10)/model.L(end,end);
%% 
% then, we take the rows associated to a few gaps (as we did in exercise A)

% Get gaps definition
r_gap = model.Input_struct.r_sens_gap;
z_gap = model.Input_struct.z_sens_gap;
t_gap = model.Input_struct.theta_sens_gap_deg;

% add some more gaps
gaps = [1 5 7 9 11 16 17 18 19 29];
rg = r_gap(gaps);
zg = z_gap(gaps);
tg = t_gap(gaps);
lg = 1; 

% Gap C matrix
i_gap = get_y_idx(model.y_type,'Gap');
i_gap = i_gap(gaps);
C_gap = model.C(i_gap,:);

% X-point C matrix
i_xp = get_y_idx(model.y_type,{'Rbound','Zbound'},1);
C_xp = model.C(i_xp,:);
%% 
% This is our choice of gaps

figure
nnodes = size(model.Input_struct.p,2);
psib = model.y_np(get_y_idx(model.y_type,'psb_c',1));
[~,hb] = plot_plasma(model.Input_struct, model.x_np(1:nnodes), psib*[1 1]);
set(hb,'linewidth',2)
hold on
xlim([0 6])
ylim([-3 4])

for i = 1 : numel(rg)
  plot([rg(i) rg(i)+lg*cosd(tg(i))], [zg(i) zg(i)+lg*sind(tg(i))],'sk-')
end
%% 
% Now, to find our transformer current pattern, we can pseudo-invert the C matrix 
% associated to the gaps and the plasma current so as to get a desired output 
% combination

dIp  = 1;
dGap = zeros(numel(rg),1);

C_It = [C_PF2Ip; 
        C_gap(:,1:10)];

Itransf = C_It\[dIp; dGap];

% And normalize
Itransf = Itransf/norm(Itransf)
%% Dynamic regulator design
% What we've seen until now is just half of the story. In fact, usually we also 
% add a dynamic regulator to our system! 
% 
% To tune our controller, we need to build our model. So, same old story; but 
% this time we'll load a pre-made model (the one we saved in the previous exercise).

% Load ss model
load PF_contr.mat PF_clsys

% Change the output equation of the PFC closed-loop system
Ip_olsys = ss(PF_clsys.a,PF_clsys.b,lm.C(i_Ip,:),lm.D(i_Ip,1:10));

% Assign I/Os
Inames = model.y_type(1:10,1);
Ip_olsys.inputname  = Inames;
Ip_olsys.outputname = 'Ip';

%% 
% This time, we'll use a PI controller. Without going into too much detail, 
% just tune the controller with the pidtune matlab command. We can tune the controller 
% on a SISO plant with the transformer current as input and the plasma current 
% as output

Ip_contr = pidtune(Ip_olsys*Itransf,'PI',3) % last parameter is the desired crossing frequency, usually the plasma current controller is not too fast
%% 
% However, we'll keep the transformer current pattern in the Ip controller for 
% later use

Ip_contr = Ip_contr*Itransf;
Ip_contr.inputname  = 'Ip';
Ip_contr.outputname = Inames
%% 
% Let's see how it works

Ip_olsys = series(Ip_contr,Ip_olsys,'name');
Ip_clsys = feedback(Ip_olsys,1);
figure
step(Ip_clsys,0:1e-3:5), grid on
save ./data/Ip_contr Ip_contr Itransf Ip_clsys