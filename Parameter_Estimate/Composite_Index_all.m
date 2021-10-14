clc
clear
close all
%% Read data form excel (No Noise) 輸入軌跡資料(將資料存在filename中的四個Excel資料夾)
filename = ["DLane2_ro.xlsx","circle2_ro.xlsx","chirp2_ro.xlsx","inv_chirp2_ro.xlsx"]; % **改**
KTrain_input_sheet = 'KTrain_input';
KTrain_output_sheet = 'KTrain_output'; 
KValidation_input_sheet = 'KValidation_input'; 
KValidation_output_sheet = 'KValidation_output'; 
bias_result_sheet = 'bias_result';
bias_tr_sheet = 'bias_tr';
ideal_tr_sheet = 'ideal_tr';

KTrain_input_data = []; KTrain_output_data = []; KValidation_input_data = []; KValidation_output_data = []; 
bias_result_data = []; bias_tr_data = []; ideal_tr_data = [];

for i = 1:length(filename)
    KTrain_input_data(:,:,i) = xlsread(filename(i),KTrain_input_sheet);  % 初始取樣訓練集輸入
    KTrain_output_data(:,:,i) = xlsread(filename(i),KTrain_output_sheet); % 初始取樣訓練集輸出
    KValidation_input_data(:,:,i) = xlsread(filename(i),KValidation_input_sheet); % 初始取樣驗證集輸入
    KValidation_output_data(:,:,i) = xlsread(filename(i),KValidation_output_sheet); % 初始取樣驗證集輸出
    bias_result_data(:,:,i) = xlsread(filename(i),bias_result_sheet);  % 具偏差參數差異量化六指標
    bias_tr_data(:,:,i) = xlsread(filename(i),bias_tr_sheet); % 具偏差參數真實系統軌跡資料
    ideal_tr_data(:,:,i) = xlsread(filename(i),ideal_tr_sheet); % 模擬模型軌跡資料
end
% for i  = 1:length(filename)
%     bias_result_data(:,:,i) = bias_result_data(:,:,i).*[1,1,1,1,0.01,0.01];
%     KTrain_output_data(:,:,i) = KTrain_output_data(:,:,i).*[1,1,1,1,0.01,0.01];
%     KValidation_output_data(:,:,i) = KValidation_output_data(:,:,i).*[1,1,1,1,0.01,0.01];
% end
% KTrain_input = KTrain_input_data;
% KTrain_output = KTrain_output_data;
% KValidation_input = KValidation_input_data;
% KValidation_output = KValidation_output_data;
% bias_result = bias_result_data;
bias_tr = bias_tr_data;
ideal_tr = ideal_tr_data;
%% 模擬模型軌跡、真實系統軌跡前處理
ideal_eq = [];
bias_eq = [];
% 等量等距插值
num = 100; % 插值數量(無須濾波，因為此部分為無雜訊軌跡)
for i = 1:length(filename)
    ideal_eq(:,:,i) = interparc(num, ideal_tr(:,1,i), ideal_tr(:,2,i),'linear');
    bias_eq(:,:,i) = interparc(num, bias_tr(:,1,i), bias_tr(:,2,i),'linear');
end
%% 原始軌跡資料繪圖
sz_r = 30;sz_s = 11;
for i = 1:length(filename)
    figure;
    hold on;
    scatter(ideal_eq(:,1,i),ideal_eq(:,2,i),sz_s,[0 0.4470 0.7410],'filled');
    scatter(bias_eq(:,1,i),bias_eq(:,2,i),sz_r,[0.8500 0.3250 0.0980],'x');
    xlabel("X-Coordinate");
    ylabel("Y-Coordinate");
%     legend
    hold off;
    axis equal;
end
%% 等量等距插值軌跡繪圖
sz_r = 30;
for i =1:length(filename)
    figure;
    hold on;
    plot(ideal_tr(:,1,i),ideal_tr(:,2,i),'color',[0 0.4470 0.7410]);
    scatter(bias_tr(:,1,i),bias_tr(:,2,i),sz_r,[0.8500 0.3250 0.0980],'x');
    % legend('Sim.','Real');
    %     title(title_name(Tr_num));
    xlabel('X-Coordinate');
    ylabel('Y-Coordinate');
    axis equal;
    hold off;
end
%% Kriging model for every index **改** 複合指標生成
bias_result=[]; KTrain_input=[]; KTrain_output=[]; KValidation_input=[]; KValidation_output=[];
lb = [48.3*0.5,50*0.5]; % upper bound **改**
ub = [48.3*1.5,50*1.5]; % lower bound **改**
fittype = 1;
SCFtype = 1;


% max min normalize 
th_sets_num = 0.005; % 篩選閥值
w_nor = [0.1,0.1,0.1,0.1,0.15,0.45]; % 各指標權重
for i = 1:length(filename)
    bias_result_n(:,:,i) = normalize(bias_result_data(:,:,i),'range');
    bias_result(:,:,i) = sum(bias_result_n(:,[1,2,3,4,5,6],i).*w_nor); % 真實系統軌跡複合指標結果

    KTrain_input = KTrain_input_data;
    KTrain_output_n(:,:,i) =  normalize(KTrain_output_data(:,:,i)','range')';
    KTrain_output(:,:,i) = sum(KTrain_output_n(:,[1,2,3,4,5,6],i).*w_nor,2); % 初始取樣訓練集複合指標結果

    KValidation_input = KValidation_input_data;
    KValidation_output_n(:,:,i) =  normalize(KValidation_output_data(:,:,i)','range')';
    KValidation_output(:,:,i) = sum(KValidation_output_n(:,[1,2,3,4,5,6],i).*w_nor,2); % 初始取樣驗證集複合指標結果
end
% 只取 Sa 與 Sc
% th_sets_num = 1;
% w = [0.7,0.3];
% for i = 1:length(filename)
% bias_result(:,:,i) = sum(bias_result_data(1,[5,6],i).*w);
% 
% KTrain_input(:,:,i) = KTrain_input_data(:,:,i);
% KTrain_output(:,:,i) = sum(KTrain_output_data(:,[5,6],i).*w,2);
% 
% KValidation_input(:,:,i) = KValidation_input_data(:,:,i);
% KValidation_output(:,:,i) = sum(KValidation_output_data(:,[5,6],i).*w,2);
% end
%% Kriging model accuracy calculation 四條軌跡各自Kriging替代模型生成；擬合結果評估(R_sq、RAAE)

% 訓練替代模形，並將訓練後的結果儲存在4Trs_Kriging.mat檔案中，已提供給具雜訊軌跡估測使用
for i = 1:length(filename)
    kparam(i) = f_variogram_fit(KTrain_input(:,:,i), KTrain_output(:,:,i), lb, ub, fittype, SCFtype);% 訓練集生成Kriging替代模形
    ykrig(:,:,i) = f_predictkrige(KValidation_input(:,:,i),kparam(i)); % 驗證集輸入，利用Kriging替代模型得到預測結果
end    
yval = KValidation_output; % 透過差異量化得到結果

% 利用ykrig、yval來評估擬合結果
N = 200;
R_sq = [];
RAAE = [];
for i = 1:length(filename) 
    R_sq(i) = 1-((sum((yval(:,:,i)-ykrig(:,:,i)).^2))/(sum((yval(:,:,i)-mean(yval(:,:,i))).^2))); 
    RAAE(i) = (sum(abs(yval(:,:,i)-ykrig(:,:,i))))/(sqrt(N*sum((yval(:,:,i)-mean(yval(:,:,i))).^2)));
end

%% contour plot **改** 各軌跡Kriging替代模型曲面繪圖
resolution = 100;
xp = [];
[xp1, xp2] = meshgrid(lb(1,1):abs(lb(1,1)-ub(1,1))/resolution:ub(1,1), ...
    lb(1,2):abs(lb(1,2)-ub(1,2))/resolution:ub(1,2));

xp(:,:,1) = xp1;
xp(:,:,2) = xp2;
xp = reshape(xp, [], 2);

fp_krig = []; fp_krig_plot = [];
for i = 1:length(filename)
    fp_krig(:,:,i) = f_predictkrige(xp,kparam(i));
    fp_krig_plot(:,:,i) = reshape(fp_krig(:,:,i), resolution+1, resolution+1);
    figure;
    surf(xp1, xp2, fp_krig_plot(:,:,i));
    xlabel("m"); % **改**
    ylabel("Iz"); % **改**
    cb = colorbar;                                     % create and label the colorbar
    cb.Label.String = 'Composite Index';
%     title_combine = strcat({'Kriging Model of '},label_name(i));
%     title(title_combine);
end

%% Parameter estimation 各軌跡替代模形估測結果
idxs_tmp = [];idxs_all = {};
th_sets = th_sets_num; % 篩選閥值
% fp_krig(:,1) = f_predictkrige(xp,kparam);
for i =1:length(filename)
    for j = 1:size(fp_krig(:,:,i),1)
        if (fp_krig(j,1,i)>=(bias_result(:,:,i)-th_sets) && fp_krig(j,1,i)<=(bias_result(:,:,i)+th_sets))
            idxs_tmp = [idxs_tmp;j];
        end
        idxs_all(:,:,i) = {idxs_tmp};
    end
    idxs_tmp = [];
end
%% 整合各軌跡有效替代模型進行參數估測，取交集得到估測結果
idxs_pass = [];
for i =1:length(filename)
    if R_sq(i) >= 0.95 && RAAE(i) <= 0.1
        idxs_pass = [idxs_pass,idxs_all{i}'];
    end
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
esti_par1 = xp(esti_idxs,1); % 取交集偏差參數m結果
esti_par2 = xp(esti_idxs,2); % 取交集偏差參數I_z結果

% 計算平均與標準差
esti_par1_mu = mean(esti_par1);
esti_par2_mu = mean(esti_par2);
esti_par1_std = std(esti_par1);
esti_par2_std = std(esti_par2);
%% 將估測結果繪圖
figure;
sz = 10;
P = [];

% 各軌跡替代模型估測結果繪圖
hold on;
P1 = scatter(xp(:,1),xp(:,2),sz,[0.7 0.7 0.7],'filled');
% c = rgb(Color,{'chocolate','coral','burlywood','cornflowerblue'});
c = [0,0.4470,0.7410;0.8500,0.3250,0.0980;0.4660,0.6740,0.1880;0.3010,0.7450,0.9330];
ct = 0;
for i = 1:length(filename)
    if R_sq(i) >= 0.95 && RAAE(i) <= 0.1
        ct = ct+1;
        P(ct) = scatter(xp(idxs_all{i},1),xp(idxs_all{i},2),sz,c(i,:),'filled');
    end
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



%% save to excel
% filename = 'composite_result_AInP(1111_15_45).xlsx';
% sheet = '4Tr';
% 
% xlswrite(filename, [esti_par1,esti_par2], sheet, 'A1');
% xlswrite(filename, [R_sq;RAAE], sheet, 'D1');
% % xlswrite(filename, krig_pass_flag, sheet, 'D6');
% % xlswrite(filename, [th_Rsq;th_RAAE],sheet,'L1');
% xlswrite(filename, [esti_par1_mu,esti_par2_mu;esti_par1_std,esti_par2_std],sheet,'D5');
% xlswrite(filename, [length(idxs_all{1}),length(idxs_all{2}),length(idxs_all{3}),length(idxs_all{4}),th_sets], sheet, 'D10');
% xlswrite(filename, [bias_result_data(:,:,1);bias_result_data(:,:,2);bias_result_data(:,:,3);bias_result_data(:,:,4)],sheet,'J1');
% xlswrite(filename, [bias_result_n(:,:,1);bias_result_n(:,:,2);bias_result_n(:,:,3);bias_result_n(:,:,4)],sheet,'J6');
% xlswrite(filename, [w_nor,bias_result(:,:,1);w_nor,bias_result(:,:,2);w_nor,bias_result(:,:,3);w_nor,bias_result(:,:,4)],sheet,'J11');