%% SOBOL INDICES ON GLOBAL SENSITIVITY FOR TADPOLE THREE WHEEL VEHICLE % **改**
clc
clear
close all
%% 輸入各軌跡的參數(參考怡平論文之最佳激發數值)
tic
total_count = 0;

% 四條軌跡名稱：軌跡參數、模擬執行時間
% DLane : 5.8253,6.5813  % 23s,  
% Circle : 2.8194,2.0000 % 11s, 
% chirp  : 0.2434,0.9905 % 30s
% inv chirp : 0.2000,0.7327 % 30s
sys_par_data = [2.8194,2.0000 ,0,1,1]; % **改**

sys_par = sys_par_data;
iii = sys_par(1);
jjj = sys_par(2);

total_count = total_count+1;
%% 1. Complex Model Sampling  **改** 利用閉迴路(open_loop)車輛模型，產生理想軌跡
% In this sampling stage, we first need to outcome the input and output
% data as the learing data for kriging fitting.
% input x: x is derived from sobol sequence depends on variable numbers and
% sampling numbers.
% output y: is the total squared error of output track and nominal values
% set.

runcycle = 100; % total sampling times % 替代模形初始取樣數量(訓練集取樣數量 + 驗證集取樣數量)
n_data_collect = 20000; % ideal trajectory num (origin 100000) % 模擬模型軌跡的資料點數量
n_real_data_collect = 100; % 真實系統軌跡的資料點數量
num = n_real_data_collect; % num 為等量等距插值數量(同真實系統軌跡資料點數量) 
tic
%---Nominal Set------------------------------------------------------
% give maneuver parameters
case_name = 'circle'; % **改** 四條軌跡切換(DLane、circle、chirp、inv_chirp)
sys_par = [iii,jjj,0,1,1]; %need a1, a2, b(better = 0), zoom_rate, zoom_length
%         sys_par = [iii,jjj];
sim_time = 11; % **改** % 四條軌跡的模擬執行時間
whole_par; % generate the nominal value of system parameters
% Drive once to generate maneuver command
[waypoint, require_velocity] = waypoints(sys_par,case_name); % create path following points
% [waypoint, require_velocity] = validation_waypoint(sys_par);
sim('tadpole_dynamic_close_loop'); % drive once, ouput its track's x and y points in cartisian coordinate 執行閉迴路架構 simulink
x_nominal = interp(linspace(0,length(x.Data),length(x.Data)),x.Data,n_data_collect);% ideal x of each points 模擬模型軌跡X資料
y_nominal = interp(linspace(0,length(y.Data),length(y.Data)),y.Data,n_data_collect);% ideal y of each points 模擬模型軌跡Y資料
% Drive command
steer_command = steer_com.Data(1,:); % 控制指令-轉向
drive_command = drive_com.Data(1,:); % 控制指令-操控
%% Ideal trajectory data process 模擬模型軌跡角度資訊序列處理、等量等距插值
ideal_tr = [x_nominal', y_nominal'];
ideal_eq = interparc(num,ideal_tr(:,1), ideal_tr(:,2),'linear');
ideal_deg360 = degree360(ideal_eq(:,1), ideal_eq(:,2));
ideal_degINRO = degreeINRO(ideal_deg360);

%% Open loop with bias m and Iz  **改** 利用開迴路架構獲得具偏差參數的真實系統軌跡
whole_par_bias; % m = 60,Iz = 60 % **改** 偏差參數調整
open_output = sim('tadpole_dynamic_open_loop1'); % 開迴路架構 simulink

x_biastr = interp(linspace(0,length(x_open.Data),length(x_open.Data)),x_open.Data,n_real_data_collect);% real x of each points 真實系統軌跡X資料
y_biastr = interp(linspace(0,length(y_open.Data),length(y_open.Data)),y_open.Data,n_real_data_collect);% real y of each points 真實系統軌跡Y資料

x_bias = x_biastr'; % 無雜訊之真實系統軌跡X資料
y_bias = y_biastr'; % 無雜訊之真實系統軌跡Y資料

% x_bias = x_biastr' + normrnd(0,0.01,[size(x_biastr',1), size(x_biastr',2)]); % 加入雜訊之真實系統軌跡X資料
% y_bias = y_biastr' + normrnd(0,0.01,[size(y_biastr',1), size(y_biastr',2)]); % 加入雜訊之真實系統軌跡Y資料

bias_tr = [x_bias, y_bias];
% equal sample
bias_eq = real_equi_smp(x_bias(:,1), y_bias(:,1), num); % 真實系統軌跡等量等距插值

th_TRAN = 0.05;th_SHAPE = 0.5; % 仿射轉換與輪廓差異量化閥值之beta_a參數
% tran sim 仿射轉換差異量化
delta_x = abs(mean(bias_eq(:,1)) - mean(ideal_eq(:,1)));
delta_y = abs(mean(bias_eq(:,2)) - mean(ideal_eq(:,2)));
lb_tran = [-delta_x*5 -delta_y*5 0 -179.9999];
ub_tran = [delta_x*5 delta_y*5 5 179.999];

fun_tran = @(x)(-TRANsim(x,ideal_eq, bias_eq(:,1),bias_eq(:,2),th_TRAN));
gs = GlobalSearch;
problem = createOptimProblem('fmincon','x0',[delta_x,delta_y,1,0],'objective',fun_tran,'lb',lb_tran,'ub',ub_tran);
[xB_tran_result,fvalB_tran_result] = run(gs,problem);

% shape sim 輪廓差異量化
bias_deg360 = degree360(bias_eq(:,1), bias_eq(:,2));
bias_degINRO = degreeINRO(bias_deg360);

fun_shape = @(x)(-SHAPEsim(x,ideal_degINRO, bias_degINRO,th_SHAPE));
rotation_value = mean(ideal_degINRO-bias_degINRO);
gs = GlobalSearch;
problem = createOptimProblem('fmincon','x0',(rotation_value),'objective',fun_shape,'lb',(-359.9999),'ub',(359.9999));
[xB_shape_result,fvalB_shape_result] = run(gs,problem);
bias_result = [xB_tran_result,fvalB_tran_result*(-1),fvalB_shape_result*(-1)]; % 真實系統軌跡之差異量化六指標
%% Plot ideal and bias trajectory 模擬模型軌跡與真實系統軌跡繪圖
figure;
hold on
plot(x_nominal,y_nominal);
sz =10;
scatter(x_bias,y_bias,sz,'filled');
xlabel('X');
ylabel('Y');
title('Ideal and Bias Trajectory');
legend('ideal','bias');
axis equal
%% Plot ideal and bias equi sample 等量等距插值後軌跡繪圖
figure;
hold on
sz =10;
scatter(ideal_eq(:,1),ideal_eq(:,2),sz);
scatter(bias_eq(:,1),bias_eq(:,2),sz,'filled');
xlabel('X');
ylabel('Y');
legend('ideal','bias');
title('Ideal and Bias Equi Sample');
axis equal
%% 在這裡利用新的參數組合去跑開迴路軌跡(open loop1) 以獲得有效 Kriging 替代模型**改**
%---Complex System Sampling------------------------------------------
% Create parameter sets from low discrepency sampling
% **改** DOE(0.5~1.5)
[m,I_x,I_y,I_z,C_alpha,C_beta,SAP,R_0,k_z,delta_fL,delta_fR,V_wind,rho,Cd,FA,mu_rolling_1,mu_rolling_2,mu_dp,mu_lp,Vx,Vy,Vz,DOE] = rand_par(runcycle);
% insert parameters into simulink model
% take the normalized number into design parameters of vehicle
fraction = 100; % each fraction has 100 samples % **改**
x_open = [];
y_open = [];
for j = 1:(runcycle/fraction)
    for i = 1:fraction
        in(i) = Simulink.SimulationInput('tadpole_dynamic_open_loop1');%drive open loop vehicle
        in(i) = in(i).setVariable('l_1',l_1);
        in(i) = in(i).setVariable('l_2',l_2);
        in(i) = in(i).setVariable('w_1',w_1);
        in(i) = in(i).setVariable('w_2',w_2);
        in(i) = in(i).setVariable('h',h);
        in(i) = in(i).setVariable('p',p);
        in(i) = in(i).setVariable('m',m(i+(j-1)*fraction));
        in(i) = in(i).setVariable('g',9.81);
        in(i) = in(i).setVariable('tire_r',0.335);
        in(i) = in(i).setVariable('tire_wdt',0.04);
        in(i) = in(i).setVariable('I_x',I_x(i+(j-1)*fraction));
        in(i) = in(i).setVariable('I_y',I_y(i+(j-1)*fraction));
        in(i) = in(i).setVariable('I_z',I_z(i+(j-1)*fraction));
        in(i) = in(i).setVariable('C_alpha',C_alpha(i+(j-1)*fraction));
        in(i) = in(i).setVariable('C_beta',C_beta(i+(j-1)*fraction));
        in(i) = in(i).setVariable('SAP',SAP(i+(j-1)*fraction));
        in(i) = in(i).setVariable('R_0',R_0(i+(j-1)*fraction));
        in(i) = in(i).setVariable('k_z',k_z(i+(j-1)*fraction));
        in(i) = in(i).setVariable('delta_fL',delta_fL(i+(j-1)*fraction));
        in(i) = in(i).setVariable('delta_fR',delta_fR(i+(j-1)*fraction));
        in(i) = in(i).setVariable('V_wind',V_wind(i+(j-1)*fraction));
        in(i) = in(i).setVariable('rho',rho(i+(j-1)*fraction));
        in(i) = in(i).setVariable('Cd',Cd(i+(j-1)*fraction));
        in(i) = in(i).setVariable('FA',FA(i+(j-1)*fraction));
        in(i) = in(i).setVariable('mu_rolling_1',mu_rolling_1(i+(j-1)*fraction));
        in(i) = in(i).setVariable('mu_rolling_2',mu_rolling_2(i+(j-1)*fraction));
        in(i) = in(i).setVariable('mu_dp',mu_dp(i+(j-1)*fraction));
        in(i) = in(i).setVariable('mu_lp',mu_lp(i+(j-1)*fraction));
        in(i) = in(i).setVariable('Vx',0);
        in(i) = in(i).setVariable('Vy',0);
        in(i) = in(i).setVariable('Vz',0);
        in(i) = in(i).setVariable('I_f_wheel',0.3);
        in(i) = in(i).setVariable('I_r_wheel',0.5);
        in(i) = in(i).setVariable('drive_command',drive_command);
        in(i) = in(i).setVariable('steer_command',steer_command);
        in(i) = in(i).setVariable('sim_time',sim_time);
    end
    % set_param('tadpole_dynamic_open_loop1','AlgebraicLoopSolver','TrustRegion')
    out = parsim(in,'ShowSimulationManager','on','ShowProgress','on');%parellel computation
    par_toc = toc;
    
    for i = 1:fraction
        x_open_pt(i,:) = interp(linspace(0,length(out(1,i).x_open.Data),length(out(1,i).x_open.Data)),out(1,i).x_open.Data,n_real_data_collect);
        y_open_pt(i,:) = interp(linspace(0,length(out(1,i).y_open.Data),length(out(1,i).y_open.Data)),out(1,i).y_open.Data,n_real_data_collect);
    end
    x_open = [x_open;x_open_pt]; % 用以擬合Kriging替代模型之初始取樣軌跡X資料
    y_open = [y_open;y_open_pt]; % 用以擬合Kriging替代模型之初始取樣軌跡Y資料
end
%% Excel data **改**
% filename = 'inv_chirp2.xlsx'; % **改**
% realX_sheet = 'real_x';
% realY_sheet = 'real_y';
% m_sheet = 'PAR_and_SIM';
% I_z_sheet = 'PAR_and_SIM';
% 
% data_x = xlsread(filename,realX_sheet);
% data_y = xlsread(filename,realY_sheet);
% data_m = xlsread(filename,m_sheet,'A1:A1000'); % **改**
% data_Iz = xlsread(filename,I_z_sheet,'B1:B1000'); % **改**
% 
% real_x = data_x;
% real_y = data_y;
% m = data_m;
% I_z = data_Iz;

%% 1000 real trajectories plot **改**
% 以下三行為 open loop 後所得到之軌跡
rng('default');
real_x = x_open'; % 用以擬合Kriging替代模型之初始取樣軌跡X資料
real_y = y_open'; % 用以擬合Kriging替代模型之初始取樣軌跡Y資料

% 加雜訊時使用，同時要去修改 sg filter
% real_x = x_open' + normrnd(0,0.01,[size(x_open',1), size(x_open',2)]);
% real_y = y_open' + normrnd(0,0.01,[size(y_open',1), size(y_open',2)]);

figure;
for i = 1:runcycle
    sz =10;
    scatter(real_x(:,i),real_y(:,i),sz,'filled');
    axis equal
    hold on
end
xlabel('X');
ylabel('Y');
title('Trajectory For Estimation');
axis equal

%% Estimation Tr's  TRAN sim 初始取樣軌跡之仿射轉換差異量化
real_eq_x = [];
real_eq_y = [];
real_eq_tmp = [];
for i = 1:size(real_x,2)
    real_eq_tmp = real_equi_smp(real_x(:,i), real_y(:,i), num);
    real_eq_x = [real_eq_x,real_eq_tmp(:,1)];
    real_eq_y = [real_eq_y,real_eq_tmp(:,2)];
    real_eq_tmp = [];
end
% 使用EU來計算兩軌跡的距離
% fval_tran_result = [];
% for i = 1:size(real_eq_x,2)
%     dis_result = DISTANTsim(ideal_eq, real_eq_x(:,i), real_eq_y(:,i));
%     fval_tran_result = [fval_tran_result;dis_result];
% end
    
th_TRAN = 0.05;
x_tran_result  = [];
fval_tran_result = [];
for i = 1:size(real_eq_x,2)
    delta_x = abs(mean(real_eq_x(:,i)) - mean(ideal_eq(:,1)));
    delta_y = abs(mean(real_eq_y(:,i)) - mean(ideal_eq(:,2)));
    % space_ls = mean(sqrt(diff(real_eq(:,1)).^2+diff(real_eq(:,2)).^2))/mean(sqrt(diff(ideal_eq(:,1)).^2+diff(ideal_eq(:,2)).^2));
    lb_tran = [-delta_x*5 -delta_y*5 0 -179.9999];
    ub_tran = [delta_x*5 delta_y*5 5 179.9999];
    
    fun_tran = @(x)(-TRANsim(x,ideal_eq, real_eq_x(:,i),real_eq_y(:,i),th_TRAN));
    gs = GlobalSearch;
    problem = createOptimProblem('fmincon','x0',[delta_x,delta_y,1,0],'objective',fun_tran,'lb',lb_tran,'ub',ub_tran);
    [x_global,fval_global] = run(gs,problem);
    x_tran_result = [x_tran_result;x_global];
    fval_tran_result = [fval_tran_result;fval_global];
end

%% Estimation Tr's SHAPE sim 初始取樣軌跡之輪廓差異量化
real_degINRO = [];
for i = 1:size(real_eq_x,2)
    real_eq_tmp = [real_eq_x(:,i),real_eq_y(:,i)];
    real_deg360_tmp = degree360(real_eq_tmp(:,1), real_eq_tmp(:,2));
    real_degINRO_tmp = degreeINRO(real_deg360_tmp);
    real_degINRO = [real_degINRO,real_degINRO_tmp];
    real_eq_tmp=[];real_deg360_tmp=[];real_degINRO_tmp=[];
end

th_SHAPE = 0.5;
x_shape_result = [];
fval_shape_result = [];
for i = 1:size(real_eq_x,2)
    fun = @(x)(-SHAPEsim(x,ideal_degINRO, real_degINRO(:,i),th_SHAPE));
    rotation_value = mean(ideal_degINRO-real_degINRO(:,i));
%     x0 = [rotation_value];
%     [x_fmincon,fval_mincon] = fmincon(fun,x0,[],[],[],[],lb_shape,ub_shape);
    gs = GlobalSearch;
    problem = createOptimProblem('fmincon','x0',(rotation_value),'objective',fun,'lb',(-359.9999),'ub',(359.9999));
    [x_global,fval_global] = run(gs,problem);
    x_shape_result = [x_shape_result;x_global];
    fval_shape_result = [fval_shape_result;fval_global];
end
%% Data collect **改**  到此部分已經得到1000個初始取樣軌跡的差異量化六個指標，接著利用此六指標
PAR_and_SIM = [];
PAR_TRAN_SIM = [];
PAR_and_SIM = [m,I_z,fval_tran_result*(-1),fval_shape_result*(-1)]; % **改**
PAR_TRAN_SIM = [m,I_z,x_tran_result,fval_tran_result*(-1),fval_shape_result*(-1)]; % **改** % 初始取樣所獲得的六個指標集合

%% Kriging model in and out **改** 建立替代模形輸入與輸出
train_num = runcycle*0.8;
validation_num = runcycle-train_num;
KTrain_input = [];
KTrain_input = [m(1:train_num),I_z(1:train_num)]; % **改**
KTrain_output = [];
KTrain_output = [x_tran_result(1:train_num,:),fval_tran_result(1:train_num)*(-1),fval_shape_result(1:train_num)*(-1)]; % **改**

KValidation_input = [];
KValidation_input = [m(train_num+1:end),I_z(train_num+1:end)]; % **改**
KValidation_output = [];
KValidation_output = [x_tran_result(train_num+1:end,:),fval_tran_result(train_num+1:end)*(-1),fval_shape_result(train_num+1:end)*(-1)]; % **改**

%% Kriging model for every index **改** 六指標各自 Kriging 替代模形建立(使用彥智學長的模型)
lb = [48.3*0.5,50*0.5]; % upper bound **改**
ub = [48.3*1.5,50*1.5]; % lower bound **改**
fittype = 1;
SCFtype = 1;

% first time fit 六個替代模形
kparam_shift_X = f_variogram_fit(KTrain_input, KTrain_output(:,1), lb, ub, fittype, SCFtype);
kparam_shift_Y = f_variogram_fit(KTrain_input, KTrain_output(:,2), lb, ub, fittype, SCFtype);
kparam_scale = f_variogram_fit(KTrain_input, KTrain_output(:,3), lb, ub, fittype, SCFtype);
kparam_rotation = f_variogram_fit(KTrain_input, KTrain_output(:,4), lb, ub, fittype, SCFtype);
kparam_TRANsim = f_variogram_fit(KTrain_input, KTrain_output(:,5), lb, ub, fittype, SCFtype);
kparam_SHAPEsim = f_variogram_fit(KTrain_input, KTrain_output(:,6), lb, ub, fittype, SCFtype);

%% Kriging model accuracy calculation 各指標替代模形擬合程度評估(R_sq、RAAE)
Kval_shift_X = f_predictkrige(KValidation_input,kparam_shift_X);
Kval_shift_Y = f_predictkrige(KValidation_input,kparam_shift_Y);
Kval_scale = f_predictkrige(KValidation_input,kparam_scale);
Kval_rotation = f_predictkrige(KValidation_input,kparam_rotation);
Kval_TRANsim = f_predictkrige(KValidation_input,kparam_TRANsim);
Kval_SHAPEsim = f_predictkrige(KValidation_input,kparam_SHAPEsim);

ykrig = [Kval_shift_X,Kval_shift_Y,Kval_scale,Kval_rotation,Kval_TRANsim,Kval_SHAPEsim];
yval = KValidation_output;
N = validation_num;
R_sq = [];
RAAE = [];
for i =1:6
    R_sq_tmp = 1-((sum((yval(:,i)-ykrig(:,i)).^2))/(sum((yval(:,i)-mean(yval(:,i))).^2))); 
    RAAE_tmp = (sum(abs(yval(:,i)-ykrig(:,i))))/(sqrt(N*sum((yval(:,i)-mean(yval(:,i))).^2)));
    R_sq = [R_sq,R_sq_tmp];
    RAAE = [RAAE,RAAE_tmp];
end

%% contour plot **改** 替代模形曲面繪圖
resolution = 100;
xp = [];
[xp1, xp2] = meshgrid(lb(1,1):abs(lb(1,1)-ub(1,1))/resolution:ub(1,1), ...
    lb(1,2):abs(lb(1,2)-ub(1,2))/resolution:ub(1,2));

xp(:,:,1) = xp1;
xp(:,:,2) = xp2;
xp = reshape(xp, [], 2);

label_name = ["Shift X","Shift Y","Scale","Rotation","TRANsim","SHAPEsim"]; % **改**
kparam = [kparam_shift_X,kparam_shift_Y,kparam_scale,kparam_rotation,kparam_TRANsim,kparam_SHAPEsim];
for i = 1:6
    fp_krig = f_predictkrige(xp,kparam(i));
    fp_krig = reshape(fp_krig, resolution+1, resolution+1);
    figure;
    surf(xp1, xp2, fp_krig);
    xlabel("m"); % **改**
    ylabel("Iz"); % **改**
    cb = colorbar;                                     % create and label the colorbar
    cb.Label.String = label_name(i);
    title_combine = strcat({'Kriging Model of '},label_name(i));
    title(title_combine);
end
%% Save to Excel **改** 資料儲存至Excel
% xlswrite('circle2_16k.xlsx',ideal_tr,'ideal_tr');
% xlswrite('circle2_16k.xlsx',bias_tr,'bias_tr');
% xlswrite('circle2_16k.xlsx',real_x,'real_x');
% xlswrite('circle2_16k.xlsx',real_y,'real_y');
% xlswrite('circle2_16k.xlsx',PAR_and_SIM,'PAR_and_SIM');
% xlswrite('circle2_16k.xlsx',PAR_TRAN_SIM,'PAR_TRAN_SIM');
% xlswrite('circle2_16k.xlsx',KTrain_input,'KTrain_input');
% xlswrite('circle2_16k.xlsx',KTrain_output,'KTrain_output');
% xlswrite('circle2_16k.xlsx',KValidation_input,'KValidation_input');
% xlswrite('circle2_16k.xlsx',KValidation_output,'KValidation_output');
% xlswrite('circle2_16k.xlsx',[R_sq;RAAE],'krig_score');
% xlswrite('circle2_16k.xlsx',bias_result,'bias_result');