%% Vertical Stabilization system
%% Preliminaries
% In this exercise, we want to design a simple VS system for the JT-60SA tokamak. 
% 
% Elongated plasmas are unstable, due to the configuration of the magnetic field 
% used to produce the elongated shape 
% 
% 
% 
% 
% 
% Hence, stabilizing the plasma is a priority in a tokamak PCS, and often a 
% dedicated actuator is used to this aim. This actuator often consists in one 
% or more _in-vessel_ circuits, which do not suffer from the slowing-down effect 
% of the tokamak walls, and thus provide a faster control action.
% 
% The actuator we'll consider here for the VS system on JT-60SA is the so-called 
% _imbalance current,_ i.e. a linear combination of the in-vessel currents (just 
% an anti-series connection) that produces an approximately radial field. Let's 
% take a look at the field lines produced by this combination of currents. We 
% can use a (linear) plasmaless response model. Let's load it

addpath ./functions
addpath ./models
addpath ./data 

plasmalessName = fullfile(pwd,'models','plasmaless_FG.mat');
plasmaless = load(plasmalessName);
%% 
% and compute the matrices of its standard ISU representation

pm.A = -plasmaless.L\plasmaless.R;
pm.B = plasmaless.L\eye(size(plasmaless.L));
pm.C = plasmaless.C;
pm.D = zeros(size(pm.C,1), size(pm.B,2)); % just a bunch of zeros
%% 
% To have a look at the field produced by the imbalance current, we choose $\delta 
% I$so that we have +1A on $I_{VSU}$ and -1A on $I_{VSL}$

dI = zeros(size(pm.C,2),1);
dI(11) = +1;
dI(12) = -1;
%% 
% then we look for the indexes of the virtual flux sensors grid (just as we 
% did in the previous exercise)

i_fg = get_y_idx(plasmaless.y_type,'Flux_grid');
n_fg = plasmaless.y_type(i_fg,1);
r_fg = plasmaless.Input_struct.r_sens(contains(plasmaless.Input_struct.names_sensors,n_fg));
z_fg = plasmaless.Input_struct.z_sens(contains(plasmaless.Input_struct.names_sensors,n_fg));
%% 
% and compute the flux variation due to the chosen $\delta I$(remember that 
% the iso-flux surfaces are also poloidal magnetic field lines!)

C_fg = pm.C(i_fg,:);
dy   = C_fg*dI;

R_fg = reshape(r_fg,30,30);
Z_fg = reshape(z_fg,30,30);
dY   = reshape(dy,30,30);

figure
plot_plasma(plasmaless.Input_struct);
hold on
contour(R_fg,Z_fg,dY,31)
xlim([ 1.2 4.5])
ylim([-3.3 3.3])
%% 
% As you can see from the picture, the field is approximately radial (and exactly 
% radial in the center of the vacuum chamber), so it is capable of providing a 
% net vertical force on the plasma ring - exactly what we need.
%% VS model
% To design our VS system, we need a model to work on. Hence, let us load our 
% linearized model (just for reference, this is a snapshot of JT-60SA Scenario 
% #2 at the Start of Flat-top)

modelName = fullfile(pwd,'models','SOF@18d66s_FG.mat'); % specify model's path
model = load(modelName); % load model variables into a structure
%% 
% and let us compute the linearized model matrices (we put them into a structure 
% for ease of reference)

lm.A = -model.L\model.R;
lm.B =  model.L\eye(size(model.L));
lm.C =  model.C;
lm.D =  zeros(size(lm.C,1), size(lm.B,2)); % again, just a bunch of zeros
lm.E = -model.L\model.LE;
lm.F =  model.F;
%% 
% Now, we need to define the inputs and outputs of our VS system. As we said, 
% for the VS system on JT-60SA we want to use the _imbalance current_ between 
% the in-vessel circuits_._ We simply choose as our control input the voltage 
% to these circuits (VSU-VSL) and connect them in anti-series.

% VSU-VSL are the 11th and 12th inputs
lm.B = lm.B(:,11)-lm.B(:,12);
lm.D = lm.D(:,11)-lm.D(:,12);
%% 
% Then, we need to choose our outputs. The outputs for the VS system are usually 
% the plasma centroid vertical velocity $\dot{z}_p(t)$ and the actuator current, 
% in our case $(I_{VSU}-I_{VSL})/2$. We start by taking the vertical position 
% and the imbalance current as follows

% Identify output positions
i_out = get_y_idx(model.y_type,{'VSU','VSL','Zpl'},1);

% Re-assign C ...
C_imbalance = (lm.C(i_out(1),:)-lm.C(i_out(2),:))/2;
C_zp        =  lm.C(i_out(3),:);
lm.C = [C_imbalance; C_zp];
% ... D
D_imbalance = (lm.D(i_out(1),:)-lm.D(i_out(2),:))/2;
D_zp        =  lm.D(i_out(3),:);
lm.D = [D_imbalance; D_zp];
% ... and F (just for consistency, we'll not use it!)
F_imbalance = (lm.F(i_out(1),:)-lm.F(i_out(2),:))/2;
F_zp        =  lm.F(i_out(3),:);
lm.F = [F_imbalance; F_zp];
%% 
% We can define a plasma state-space object as

VS_sys = ss(lm.A,lm.B,lm.C,lm.D,...
            'inputname',{'VSU-VSL'},'outputname',{'(IVSU-IVSL)/2','Zp'});
%% 
% and verify that our model has an unstable eigenvalue (aka a pole in the right 
% half of the complex plane)

max(real(eig(VS_sys)))
%% 
% Incidentally, the next few eigenvalues have 0 real part, and are due to the 
% superconducting coils (which theoretically behave as integrators)

esort(real(eig(VS_sys)))
%% 
% We'll carry out the design of the controller on a simple case, where we consider 
% an ideal derivative of the vertical plasma position. Since, in our model, we 
% have
% 
% $$\dot{x}(t)=Ax(t)+Bu(t) \\y(t) = Cx(t)+Du(t)$$
% 
% we can compute the derivative of y as
% 
% $$\dot{y}(t) = CAx(t)+CBu(t) + D\dot{u}(t)$$
% 
% In our case, D is just a bunch of zeros, so we can neglect the last term. 
% Applying this procedure to the second output of our model, we have

VS_olsys = ss(VS_sys.a, ...
            VS_sys.b, ...
            [VS_sys.c; VS_sys.c(2,:)*VS_sys.a], ...
            [VS_sys.d; VS_sys.c(2,:)*VS_sys.b]);
%% 
% Again, assign I/O names (this extra step will turn out to be quite handy in 
% a moment)

VS_olsys.outputname = {'(IVSU-IVSL)/2','Zp','Zdot'};
VS_olsys.inputname  = {'VSU-VSL'};
%% VS design
% Next step is to define a test case for our constrol system. A worst-case scenario 
% often use to validate a VS design is to consider a so-called VDE (Vertical Displacement 
% Event), i.e. a displacement along the unstable current pattern which results 
% into a given $\delta z_p$. 
% 
% First we need to get our unstable direction, which is the eigenvector of the 
% A matrix associated to the unstable eigenvalue

[V,D] = eig(lm.A);
[gamma,i_gamma] = max(real(diag(D))); % gamma is the growth-rate of our unstable mode
x_vde = V(:,i_gamma); % normalized to 1A
%% 
% The (very small) plasma vertical displacement associated to this unstable 
% pattern is

z0 = lm.C(2,:)*x_vde
%% 
% so, if we ask (for instance) a 2cm vde

x0 = x_vde*2e-2/z0;
figure
subplot(221)
bar(x0(1:10))
title('PF currents')
subplot(222)
bar(x0(11:12))
title('IVSU-IVSL')
subplot(223)
bar(x0(13:end-1))
title('passive currents')
subplot(224)
bar(x0(end))
title('Ip')
%% 
% There are different ways to design a VS system. Let's go for a very simple 
% one. First, of all, we fix the structure of our controller as
% 
% $$u = K_1 I_{VS} + K_2 \dot{z}_p(t)$$
% 
% (notice that in general the velocity gain is scaled by the scenario plasma 
% current, $K_2 = \bar{K}_2 I_{p0}$).
% 
% We will then use the optimization tools provided by matlab. The objective 
% function is defined in the file VS_objf.m, and is just a weighted sum of the 
% squares of the imbalance current and the plasma vertical velocity samples obtained 
% from a 0.5s simulation. We look for its minimum

Ts    = 1e-3; % sampling time for the control system
Tspan = (0:Ts:0.5)';

% Just some options for the solver
options = optimset('MaxIter',1000,'MaxFunEvals',1500,'Display','Iter');

% Solve the optimization problem
[K_opt,J_opt] = fmincon(@(k)VS_objf(k,VS_olsys,x0,Tspan),[-0.03 800],[],[],[],[],[],[],@(k)VS_nlcon(k,VS_olsys),options);
%% 
% Finally, we define our controller as a (static) ss object. Notice that we 
% can use the I/O names to automatically make the connections between the blocks!

VS_contr = ss(K_opt,'inputname',{'(IVSU-IVSL)/2','Zdot'},'outputname',{'VSU-VSL'});
%% Test
% Test #1
% Now it's time to test our controller. Define the closed-loop system by connecting 
% the controller to the open-loop model of our plant

VS_clsys = feedback(VS_olsys,VS_contr,'name'); % the 'names' option indicates that the connections are made automatically based on the I/O tags
%% 
% Then simulate and plot some figures

[y,t] = lsim(VS_clsys,Tspan*0,Tspan,x0,'zoh');

figure
subplot(311)
plot(t,y(:,1))
title(VS_olsys.outputname{1})
grid on
subplot(312)
plot(t,y(:,2))
title(VS_olsys.outputname{2})
grid on
subplot(313)
plot(t,y(:,3))
title(VS_olsys.outputname{3})
grid on
%% 
% Good: the velocity is brought to 0 very rapidly, which is what we wanted. 
% Perhaps even a bit _too_ rapidly - we'll discuss this in a moment.
% Test #2 (a cautionary tale)
% In general, not only we need a _fast_ controller, but also a _robust_ one. 
% Control theory provides different ways to check the robustness properties of 
% our closed-loop system, as we want to make sure that our controller works even 
% when disturbances/uncertainties come at play. In this case, we can take a look 
% at the Nyquist plot of our solution

figure
nyquist(VS_contr*VS_olsys([1 3],:))
%% 
% Ugly! Let's take a closer look

figure
h=nyquistplot(VS_contr*VS_olsys([1 3],:));
opt = getoptions(h);
opt.XLim = {[-30, 10]};
opt.YLim = {[-30, 30]};
setoptions(h,opt);
%% 
% We are _extremely_ close to instability! 
% 
% In general, the VS system is a quite fast/delicate controller, for which the 
% power-supplies could have a significant impact. For this reason, it's a good 
% idea to take into account the power supplies of our in-vessel coils explicitly. 
% The model we'll consider for JT-60SA's power supplies is the following (1st 
% order dynamics + pure delay)

PS_sys = tf(1,[3e-3,1]);    % 1st order dynamics
PS_sys.inputDelay = 1.5e-3; % pure delay
PS_sys.inputname  = {'VSU-VSL'};
PS_sys.outputname = {'VSU-VSL'}; % just use the same name
%% 
% Moreover, usually the derivative of the vertical position of the plasma is 
% obtained numerically by means of a High-Pass Filter (HPF)

HPF = tf([1,0], [1e-3,1])
%% 
% We can verify that this is indeed a derivative filter, at least for (relatively) 
% low frequency signals. Consider the following example

x = 0:1e-2:5;
y = sin(2*x);
% y = sin(1000*x); % This case will NOT work (the frequency is too high!)

dy = lsim(HPF,y,x);
figure
plot(x,y,'b',x,dy,'x-r',x(1:end-1), diff(y)./diff(x), '.k')
legend({'signal','filter','numerical derivative'},'location','best')
%% 
% Moreover, for the simulation we'll need to feed-through the IV current, so 
% let's turn our HPF into a matrix

HPF = [1 0; ...
       0 1; ...
       0 HPF]; 
HPF.inputname  = {'(IVSU-IVSL)/2','Zp'};
HPF.outputname = {'(IVSU-IVSL)/2','Zp','Zdot'};
%% 
% We can now connect our systems (OL just stands for 'open-loop' here)

VS_olsys2 = series(PS_sys, VS_sys,'name');
VS_olsys2 = series(VS_olsys2,HPF,   'name'); 
%% 
% And see what happens now with the controller we had found before

VS_clsys2 = feedback(VS_olsys2,VS_contr,'name');

VS_clsys2.InputDelay = VS_clsys2.InternalDelay;
VS_clsys2.InternalDelay = 0;
[y,t] = lsim(VS_clsys2,Tspan*0,Tspan,[0;x0;0],'zoh');

figure
subplot(311)
plot(t,y(:,1))
title(VS_olsys.outputname{1})
grid on
subplot(312)
plot(t,y(:,2))
title(VS_olsys.outputname{2})
grid on
subplot(313)
plot(t,y(:,3))
title(VS_olsys.outputname{3}) 
grid on
%% 
% boom! It broke...
% 
% Sometimes, very good performance means very poor robustness :(
% 
% Let's look again at our Nyquist plot

figure
h=nyquistplot(VS_contr*VS_olsys2([1 3],:));
opt = getoptions(h);
opt.XLim = {[-30, 10]};
opt.YLim = {[-30, 30]};
setoptions(h,opt);
%% 
% As you can see, now we have one encirclement in the cw direction. The "spirals" 
% close to the origin are due to the delay introduced by the power supplies. 
% 
% We can escape this circle by reducing the controller gain. Let's try a "quieter" 
% controller

VS_contr2 = VS_contr*0.05; % much slower!
VS_clsys2 = feedback(VS_olsys2,VS_contr2,'name');

figure
h=nyquistplot(VS_contr2*VS_olsys2([1 3],:));
opt = getoptions(h);
opt.XLim = {[-5, 5]};
opt.YLim = {[-10, 10]};
setoptions(h,opt);
%% 
% We escaped the smaller loop and now the controller should work. Let's try 
% it in simulation

VS_clsys2.InputDelay = VS_clsys2.InternalDelay;
VS_clsys2.InternalDelay = 0; % this is just a quick'n'dirty trick to make lsim work
[y,t] = lsim(VS_clsys2,Tspan*0,Tspan,[0;x0;0],'zoh');

figure
subplot(311)
plot(t,y(:,1))
grid on
title(VS_olsys.outputname{1})

subplot(312)
plot(t,y(:,2))
grid on
title(VS_olsys.outputname{2}) 

subplot(313)
plot(t,y(:,3))
grid on 
title(VS_olsys.outputname{3}) % ...but safer :)

%% 
% Not bad. Finally, we save the controller for later use

VS_contr = VS_contr2;
save ./data/VS_contr VS_contr
%% 
%