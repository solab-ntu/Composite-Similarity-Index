%% Random parameters
function [m,I_x,I_y,I_z,C_alpha,C_beta,SAP,R_0,k_z,delta_fL,delta_fR,V_wind,rho,Cd,FA,mu_rolling_1,mu_rolling_2,mu_dp,mu_lp,Vx,Vy,Vz,DOE_out] = rand_par(runcycle)
%% Sobol Sequence Sampling
Q = sobolset(2,'skip',200000);% create a different sobol sequence
% DOE = net(Q,runcycle);
DOE = ((net(Q,runcycle)*2-1)*0.5)+1;
DOE_out = DOE;
%% unknown parameters
% m = 48.3*ones(runcycle,1);
m = 48.3*DOE(:,1);
I_x = 12.1*ones(runcycle,1);
I_y = 36.7*ones(runcycle,1);
I_z = 50*DOE(:,2);
% I_z = 50*ones(runcycle,1);
% I_z = 60*ones(runcycle,1);

% tire
C_alpha = 1000*ones(runcycle,1);
% C_alpha = 1000*DOE(:,2);
% C_beta = 50*DOE(:,3);
C_beta = 50*ones(runcycle,1);
% SAP = -0.015+0.03*DOE(:,4);
SAP = 0.01*ones(runcycle,1);
R_0 = 0.335*ones(runcycle,1);
k_z = 30000*ones(runcycle,1);

%% uncertainties
% steer
delta_fL = -1*ones(runcycle,1);
delta_fR = 3.5*ones(runcycle,1);

%% resistance
% aero
V_wind = 0*ones(runcycle,1);
rho = 1.225*ones(runcycle,1);
Cd = 0.25*ones(runcycle,1);
FA = 1*ones(runcycle,1);

%rolling resistance
% mu_rolling_1 = 0.009*DOE(:,5);
mu_rolling_1 = 0.009*ones(runcycle,1);
% mu_rolling_2 = 0.0005*DOE(:,6);
mu_rolling_2 = 0.0005*ones(runcycle,1);
mu_dp = 0.8*ones(runcycle,1);
mu_lp = 0.8*ones(runcycle,1);

%% initial value
Vx = 0;
Vy = 0;
Vz = 0;





