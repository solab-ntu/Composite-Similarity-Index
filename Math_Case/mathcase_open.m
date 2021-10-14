% 論文數學案例(開放軌跡)
clc;
clear;
close all;

%% 數學範例軌跡(原始sin、位移旋轉之sin、振幅改變之sin)
% sim trajectory : y=sin(x) 理想軌跡
x_s = linspace(0,6*pi,1000);
y_s = sin(x_s);
sim_tr = [x_s',y_s'];

% real trajectory 1 : 旋轉 = 5;位移=(3,3) 真實軌跡1
x = linspace(0,6*pi,50);
y = sin(x);
z = [3,3,1,5];
ori_tr = [x',y'];
angle = deg2rad(z(4));
x_r1 = z(3)*(cos(angle)*ori_tr(:,1) - sin(angle)*ori_tr(:,2)) + z(1);
y_r1 = z(3)*(sin(angle)*ori_tr(:,1) + cos(angle)*ori_tr(:,2)) + z(2);
real_tr1 = [x_r1, y_r1];

% real trajectory 1 : y=2sin(x) 振幅改變真實軌跡2
x_r2 = linspace(0,6*pi,50);
y_r2 = 2*sin(x_r2);
real_tr2 = [x_r2', y_r2'];

% trajectory plot
sz =15;
sz_s = 30;

figure
hold on
plot(sim_tr(:,1), sim_tr(:,2),'-','color','k');%[0.8500 0.3250 0.0980]
plot(real_tr1(:,1), real_tr1(:,2),'x','color',[0 0.4470 0.7410]);%[0 0.4470 0.7410]
plot(real_tr2(:,1), real_tr2(:,2),'s','color',[0.8500 0.3250 0.0980]);%[0 0.4470 0.7410]
scatter(sim_tr(1,1), sim_tr(1,2),sz_s,'k','filled');
scatter(real_tr1(1,1), real_tr1(1,2),sz,[0 0.4470 0.7410],'filled');
scatter(real_tr2(1,1), real_tr2(1,2),sz,[0.8500 0.3250 0.0980],'filled');
xlabel('X-Coordinate');
ylabel('Y-Coordinate');
legend('$T^s$','$T^{r}_1$','$T^{r}_2$','interpreter','latex');
% annotation('arrow', [1.3 .5], [.6 .5]);
grid on
hold off
axis equal;
%% real trajectory 前處理
% 理想軌跡前處理
sim_eq = interparc(50, sim_tr(:,1), sim_tr(:,2),'linear'); % 等量等距插值
sim_deg360 = degree360(sim_eq(:,1), sim_eq(:,2)); % 計算tan角度
sim_degINRO = degreeINRO(sim_deg360); % 角度序列再處理

% 真實系統軌跡1前處理
real1_eq = interparc(50,real_tr1(:,1), real_tr1(:,2),'linear'); % 等量等距插值
real1_deg360 = degree360(real1_eq(:,1), real1_eq(:,2)); % 計算tan角度
real1_degINRO = degreeINRO(real1_deg360); % 角度序列再處理

% 真實系統軌跡2前處理
real2_eq = interparc(50,real_tr2(:,1), real_tr2(:,2),'linear'); % 等量等距插值
real2_deg360 = degree360(real2_eq(:,1), real2_eq(:,2)); % 計算tan角度
real2_degINRO = degreeINRO(real2_deg360); % 角度序列再處理
 
% trajectory plot
sz =15;
sz_s = 30;
figure
hold on
plot(sim_eq(:,1), sim_eq(:,2),'.','color','k');
plot(real1_eq(:,1), real1_eq(:,2),'x','color',[0 0.4470 0.7410]);
plot(real2_eq(:,1), real2_eq(:,2),'s','color',[0.8500 0.3250 0.0980]);
scatter(sim_eq(1,1), sim_eq(1,2),sz_s,'k','filled');
scatter(real1_eq(1,1), real1_eq(1,2),sz,[0 0.4470 0.7410],'filled');
scatter(real2_eq(1,1), real2_eq(1,2),sz,[0.8500 0.3250 0.0980],'filled');
xlabel('X-Coordinate');
ylabel('Y-Coordinate');
legend('$T^s$','$T^{r}_1$','$T^{r}_2$','interpreter','latex');
% annotation('arrow', [1.3 .5], [.6 .5]);
grid on
hold off
axis equal;
%%  real trajectory 1 真實系統軌跡1，仿射與輪廓差異量化(最佳化)
delta_x1 = abs(mean(real1_eq(:,1)) - mean(sim_eq(:,1)));
delta_y1 = abs(mean(real1_eq(:,2)) - mean(sim_eq(:,2)));
lb_tran1 = [-delta_x1*5 -delta_x1*5 0 -5];
ub_tran1 = [delta_y1*5 delta_y1*5 5 5];
% lb_tran1 = [-5 -5 0 -5];
% ub_tran1 = [5 5 5 5];

% Affine transformation similarity 仿射轉換差異量化
fun_tran = @(x)(-TRANsim(x,sim_eq, real1_eq(:,1),real1_eq(:,2),0.05));
gs = GlobalSearch;
problem = createOptimProblem('fmincon','x0',[0,0,1,0],'objective',fun_tran,'lb',lb_tran1,'ub',ub_tran1);
[xB_tran_result1,fvalB_tran_result1] = run(gs,problem);

% Shape similarity 輪廓差異量化
fun_shape = @(x)(-SHAPEsim(x,sim_degINRO,real1_degINRO,1));
rotation_value = mean(sim_degINRO-real1_degINRO);
gs = GlobalSearch;
problem = createOptimProblem('fmincon','x0',(rotation_value),'objective',fun_shape,'lb',(-359.999),'ub',(359.9999));
[x_global1,fval_global1] = run(gs,problem);

bias_result_1 = round([xB_tran_result1,fvalB_tran_result1*(-1),fval_global1*(-1)],3);
w_nor = [0.1,0.1,0.1,0.1,0.15,0.45];
bias_result_n1 = normalize(bias_result_1,'range');
cmp_index1 = sum(bias_result_n1([1,2,3,4,5,6]).*w_nor);

%%  real trajectory 2 真實系統軌跡2，仿射與輪廓差異量化(最佳化)
delta_x2 = abs(mean(real2_eq(:,1)) - mean(sim_eq(:,1)));
delta_y2 = abs(mean(real2_eq(:,2)) - mean(sim_eq(:,2)));
lb_tran2 = [-delta_x2*5 -delta_x2*5 0 -5];
ub_tran2 = [delta_y2*5 delta_y2*5 5 5];
% lb_tran2 = [-5 -5 0 -5];
% ub_tran2 = [5 5 5 5];

% Affine transformation similarity 仿射轉換差異量化
fun_tran = @(y)(-TRANsim(y,sim_eq, real2_eq(:,1),real2_eq(:,2),0.05));
gs = GlobalSearch;
problem = createOptimProblem('fmincon','x0',[0,0,1,0],'objective',fun_tran,'lb',lb_tran2,'ub',ub_tran2);
[xB_tran_result2,fvalB_tran_result2] = run(gs,problem);

% Shape similarity 輪廓差異量化
fun_shape = @(x)(-SHAPEsim(x,sim_degINRO,real2_degINRO,1));
rotation_value = mean(sim_degINRO-real2_degINRO);
gs = GlobalSearch;
problem = createOptimProblem('fmincon','x0',(rotation_value),'objective',fun_shape,'lb',(-359.999),'ub',(359.9999));
[x_global2,fval_global2] = run(gs,problem);

bias_result_2 = round([xB_tran_result2,fvalB_tran_result2*(-1),fval_global2*(-1)],3);
w_nor = [0.1,0.1,0.1,0.1,0.15,0.45];
bias_result_n2 = normalize(bias_result_2,'range');
cmp_index2 = sum(bias_result_n2([1,2,3,4,5,6]).*w_nor);

%% 轉換後結果繪圖
z = xB_tran_result2;
angle = deg2rad(z(4));
x_tran = z(3)*(cos(angle)*sim_eq(:,1) - sin(angle)*sim_eq(:,2)) + z(1);
y_tran = z(3)*(sin(angle)*sim_eq(:,1) + cos(angle)*sim_eq(:,2)) + z(2);
figure
hold on
plot(x_tran, y_tran);
plot(real2_eq(:,1), real2_eq(:,2),':x','color',[0 0.4470 0.7410]);
plot(sim_eq(:,1), sim_eq(:,2),':s','color',[0.8500 0.3250 0.0980]);
axis equal
hold off
