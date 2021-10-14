%% fixed parameters
% chassis
l_1 = 0.724;
l_2 = 0.719;
w_1 = 0.3065;
w_2 = 0.3135;
h = 0.4896;
p = 1.43;
% mass and  inertia
m = 60; % 48.3
g = 9.81;
tire_r = 0.335;
tire_wdt = 0.04; %tire wide /2

%% unknown parameters
I_x = 12.1;
I_y = 36.7;
I_z = 60; % 50

I_f_wheel = 0.3;
I_r_wheel = 0.5;

% tire


% C_alpha = 200;
% C_beta = 0.04;
C_alpha = 1000; %1000
C_beta = 50;
SAP = 0.01;
R_0 = 0.335;
k_z = 30000;



%% uncertainties
% steer
delta_fL = -1;
delta_fR = 3.5;

%% resistance
% aero
V_wind = 0;
rho = 1.225;
Cd = 0.25;
FA = 1;

%rolling resistance
mu_rolling_1 = 0.009;
mu_rolling_2 = 0.0005;
mu_dp = 0.8;
mu_lp = 0.8;

%% initial value
Vx = 0;
Vy = 0;
Vz = 0;