clc
clear
close all
%% 將四軌跡替代模形儲存在 4Trs_Kriging.mat 中
pass_krig = [1,3,4]; % 獲得通過驗證的軌跡替代模形
k = load('4Trs_Kriging.mat');
kparam = k.kparam(pass_krig);
%% load datas 輸入無雜訊軌跡資料(將資料存在 filename 中的四個Excel資料夾)
filename_all = ["DLane2_ro.xlsx","circle2_ro.xlsx","chirp2_ro.xlsx","inv_chirp2_ro.xlsx"]; % **改**
filename = filename_all(pass_krig);
ideal_tr_sheet = 'ideal_tr';


ideal_tr_data = [];
for i =1:length(filename)
    ideal_tr_data(:,:,i) = xlsread(filename(i),ideal_tr_sheet);
end
ideal_tr = ideal_tr_data;

% 加入具雜訊的軌跡(雜訊大小為sigma=0.1，且具雜訊的四條操作軌跡儲存在 Tr4_noise0.1 Excel 檔案中)
rng('default');
x_bias = []; y_bias=[];
x_bias_data =  xlsread("Tr4_noise0.1.xlsx","noise_tr_x") 

y_bias_data =  xlsread("Tr4_noise0.1.xlsx","noise_tr_y") ;
x_bias = x_bias_data(:,pass_krig);
y_bias = y_bias_data(:,pass_krig);
%% 原始軌跡資料繪圖
sz_r = 30; title_name = ["DLane Trajectory","Circle Trajectory","Chirp Trajectory","Inv Chirp Trajectory"];

for i = 1:length(filename)
    figure;
    hold on;
    p_i = plot(ideal_tr(:,1,i),ideal_tr(:,2,i),'color',[0 0.4470 0.7410]);
    s_b = scatter(x_bias(:,i),y_bias(:,i),sz_r,[0.8500 0.3250 0.0980],'x');
    legend('Sim.','Real');
    %     title(title_name(Tr_num));
    xlabel('X-Coordinate');
    ylabel('Y-Coordinate');
    axis equal;
    hold off;
end
%%  模擬模型軌跡、真實系統軌跡
num = 100 ; % 等量等距差值數量
ideal_eq = []; ideal_deg360 = []; ideal_degINRO = [];

% 模擬模型軌跡等量等距插值
for i = 1:length(filename)
    ideal_eq(:,:,i) = interparc(num,ideal_tr(:,1,i), ideal_tr(:,2,i),'linear');
    ideal_deg360(:,i) = degree360(ideal_eq(:,1,i), ideal_eq(:,2,i));
    ideal_degINRO(:,i) = degreeINRO(ideal_deg360(:,i));
end
bias_eq = [];
bias_eq_plot = [];

th_com = [];
th_TRAN = [0.14, 0.14, 0.14]; % 調整出之仿射差異量化閥值參數 beta_a
th_SHAPE = [4.85,3.1, 3.6]; % 調整出之輪廓差異量化閥值參數 beta_a
[th_com1, th_com2] = meshgrid(th_TRAN,th_SHAPE);

th_com(:,:,1) = th_com1;
th_com(:,:,2) = th_com2;
th_com = reshape(th_com, [], 2);
%% 真實系統軌跡差異量化
for i = 1:length(filename)
    bias_eq = real_equi_smp(x_bias(:,i), y_bias(:,i), num); % 濾波後等量等距插值
    bias_eq_plot(:,:,i) = bias_eq;
    bias_deg360 = degree360(bias_eq(:,1), bias_eq(:,2));
    bias_degINRO = degreeINRO(bias_deg360);
    
    delta_x = abs(mean(bias_eq(:,1)) - mean(ideal_eq(:,1,i)));
    delta_y = abs(mean(bias_eq(:,2)) - mean(ideal_eq(:,2,i)));
    lb_tran = [-delta_x*5 -delta_y*5 0 -179.9999];
    ub_tran = [delta_x*5 delta_y*5 5 179.999];


    % TRANsim 仿射轉換差異量化
    fun_tran = @(x)(-TRANsim(x,ideal_eq(:,:,i), bias_eq(:,1),bias_eq(:,2),th_TRAN(i)));
    gs = GlobalSearch;
    problem = createOptimProblem('fmincon','x0',[delta_x,delta_y,1,0],'objective',fun_tran,'lb',lb_tran,'ub',ub_tran);
    [xB_tran_result,fvalB_tran_result] = run(gs,problem);

    % SHAPEsim 輪廓差異量化
    fun_shape = @(x)(-SHAPEsim(x,ideal_degINRO(:,i), bias_degINRO,th_SHAPE(i)));
    rotation_value = mean(ideal_degINRO(:,i)-bias_degINRO);
    gs = GlobalSearch;
    problem = createOptimProblem('fmincon','x0',(rotation_value),'objective',fun_shape,'lb',(-359.9999),'ub',(359.9999));
    [xB_shape_result,fvalB_shape_result] = run(gs,problem);
    bias_result_data(i,:) = [xB_tran_result,fvalB_tran_result*(-1),fvalB_shape_result*(-1)];
end
%% 等量等距插值軌跡繪圖
for i = 1:length(filename)
    figure;
    sz_r  = 12;
    hold on
    scatter(ideal_eq(:,1,i),ideal_eq(:,2,i),sz_r,[0 0.4470 0.7410],'filled');
    scatter(bias_eq_plot(:,1,i),bias_eq_plot(:,2,i),sz_r,[0.8500 0.3250 0.0980],'x');
    xlabel('X-Coordinate');
    ylabel('Y-Coordinate');
    legend('Sim.','Real');
    axis equal;
    hold off
end
%% 複合指標計算並儲存在bias_result中
lb = [48.3*0.5,50*0.5]; % upper bound **改**
ub = [48.3*1.5,50*1.5]; % lower bound **改**
fittype = 1;
SCFtype = 1;

% max min normalize 
th_sets_num = 0.005;
w_nor = [0.1,0.1,0.1,0.1,0.15,0.45]; % 權重設計
bias_result_n = []; bias_result=[];
for i = 1:length(filename)
    bias_result_n(i,:) = normalize(bias_result_data(i,:),'range');
    bias_result(i) = sum(bias_result_n(i,[1,2,3,4,5,6]).*w_nor);
end
% 只取 Sa 與 Sc
% th_sets_num = 3;
% w = [0.3,0.7]; % [1,1,1,1,3,3]
% for i = 1:length(filename)
%     bias_result(i) = sum(bias_result_data(i,[5,6]).*w);
% end
%% 使用無雜訊軌跡有效之替代模型進行參數估測(各自軌跡估測結果)
resolution = 100;
xp = [];
[xp1, xp2] = meshgrid(lb(1,1):abs(lb(1,1)-ub(1,1))/resolution:ub(1,1), ...
    lb(1,2):abs(lb(1,2)-ub(1,2))/resolution:ub(1,2));

xp(:,:,1) = xp1;
xp(:,:,2) = xp2;
xp = reshape(xp, [], 2);


idxs_tmp = [];fp_krig = [];idxs_all={};
th_sets = th_sets_num; % 估測閥值設計
for i = 1:length(filename)
    fp_krig(:,i) = f_predictkrige(xp,kparam(i));
end
for i = 1:length(filename)
    for j = 1:size(fp_krig(:,i),1)
        if (fp_krig(j,i)>=(bias_result(i)-th_sets) && fp_krig(j,i)<=(bias_result(i)+th_sets))
            idxs_tmp = [idxs_tmp;j];
        end
    end
    idxs_all(:,:,i) = {idxs_tmp};
    idxs_tmp = [];
end
%% 整合各軌跡有效替代模型進行參數估測，取交集得到估測結果
idxs_pass = [];
for i =1:length(filename)
    idxs_pass = [idxs_pass,idxs_all{i}'];
end
bin_edge = unique(idxs_pass);
count_arr = histc(idxs_pass,bin_edge);
count_st = max(count_arr);
% count_nd = max(count_arr(count_arr<count_st));
% count_rd = max(count_arr(count_arr<count_nd));

index_st = find(count_arr == count_st);
% index_nd = find(count_arr == count_nd);
% index_rd = find(count_arr == count_rd);

esti_idxs = bin_edge(index_st);
esti_par1 = xp(esti_idxs,1);
esti_par2 = xp(esti_idxs,2);
esti_par1_mu = mean(esti_par1);
esti_par2_mu = mean(esti_par2);
esti_par1_std = std(esti_par1);
esti_par2_std = std(esti_par2);
%% 將估測結果繪圖

% 各軌跡替代模型估測結果繪圖
figure;
sz = 10;
P = [];
hold on;
P1 = scatter(xp(:,1),xp(:,2),sz,[0.7 0.7 0.7],'filled');
% c = rgb(Color,{'chocolate','coral','burlywood','cornflowerblue'});
c = [0,0.4470,0.7410;0.4660,0.6740,0.1880;0.3010,0.7450,0.9330];
for i = 1:length(filename)
    P(i) = scatter(xp(idxs_all{i},1),xp(idxs_all{i},2),sz,c(i,:),'filled');
end
P2 = scatter(xp(esti_idxs,1),xp(esti_idxs,2),sz,[1 1 0],'filled');
sz_x = 50;
P3 = scatter(60,60,sz_x,[1,0,0],'x','linewidth',1.3);
label_name = ["DLane","Chirp","Inv Chirp"]; %"Circle"
plgend = label_name;
p_all = [P2,P,P1,P3]; 
lgd = legend(p_all,['Estimation Points',plgend,'Initial Space','Bias Value'],'location','eastoutside','Interpreter','latex');
lgd.FontSize = 8;
xlabel("m",'Interpreter','latex'); % **改**
ylabel("Iz",'Interpreter','latex'); % **改**
axis equal;
hold off

% 整合各軌跡有效替代模形估測結果繪圖
for i = 1:length(filename)
    figure;
    hold on;
    p1 = scatter(xp(:,1),xp(:,2),sz,[0.7 0.7 0.7],'filled');
    p2 = scatter(xp(idxs_all{i},1),xp(idxs_all{i},2),sz,[1 1 0],'filled');
    sz_x = 50;
    p3 = scatter(60,60,sz_x,[1,0,0],'x','linewidth',1.3);
    p_all = [p2,p1,p3];
    lgd = legend(p_all,'Estimation Points','Initial Space','Bias Value','location','eastoutside','Interpreter','latex');
    lgd.FontSize = 8;
    xlabel("m",'Interpreter','latex'); % **改**
    ylabel("Iz",'Interpreter','latex'); % **改**
    axis equal;
    hold off;
end
%% 儲存資料
% filename = 'composite_result_AInP(1111_15_45).xlsx';
% sheet = 'no_4Tr';
% % 
% xlswrite(filename, [esti_par1,esti_par2], sheet, 'A1');
% % xlswrite(filename, [R_sq;RAAE], sheet, 'D1');
% % xlswrite(filename, krig_pass_flag, sheet, 'D6');
% % xlswrite(filename, [th_Rsq;th_RAAE],sheet,'L1');
% xlswrite(filename, [esti_par1_mu,esti_par2_mu;esti_par1_std,esti_par2_std],sheet,'D5');
% xlswrite(filename, [length(idxs_all{1}),length(idxs_all{2}),length(idxs_all{3}),th_sets], sheet, 'D10');
% xlswrite(filename, [th_TRAN',th_SHAPE'], sheet, 'D13');
% xlswrite(filename, [bias_result_data(1,:);bias_result_data(2,:);bias_result_data(3,:)],sheet,'J1');
% xlswrite(filename, [bias_result_n(1,:);bias_result_n(2,:);bias_result_n(3,:)],sheet,'J6');
% xlswrite(filename, [w_nor,bias_result(1);w_nor,bias_result(2);w_nor,bias_result(3)],sheet,'J11');