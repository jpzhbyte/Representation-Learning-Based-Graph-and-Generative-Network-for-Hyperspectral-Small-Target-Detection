close all;
%% 数据
load('.\Data\HYDICE_data.mat');
load('.\Data\initial_det.mat');
%% 加载feature
[row,col,bands]=size(data);
A = ones(row,col);
load('.\Features\gcnfeature.mat');
features = hyperNormalize(features);
TE2d = hyperConvert2d(TE);
x = find(TE2d>0);
Pred_TE = zeros(size(TE2d));
Pred_TE = max(features', [], 1);
Pred_TE3d = hyperConvert3d(Pred_TE,row,col);
Pred_CM = Pred_TE3d;
Pred_CM = hyperNormalize(Pred_CM);
%% 测试网络输出
num_set = [1:199];
[r,c,b] = size(data_h2);
a_Result = zeros(200,2);
for i = 1:15
    test = load(strcat('.\test\HYDICE\HYDICE_',num2str(num_set(i)),'.mat'));
    t_0 = struct2array(test);
    t = hyperNormalize(t_0);
    t = double(t);
    t_1=reshape(t,col,row,bands);
    t_3=permute(t_1,[2,1,3]);
    X = reshape(t_3,col*row,bands);
    X = X';
    [N M] = size(X);
    X_mean = mean(X,2);
    X = X - repmat(X_mean, [1 M]);
    Sigma = (X * X')/(M-1);
    Sigma_inv = inv(Sigma);
    D=zeros(1,M);
    for m = 1:M
        D(m) = X(:, m)' * Sigma_inv * X(:, m);
    end
    result_RX1 = reshape(D,col,row);
    result_RX2 = result_RX1';
    result = hyperNormalize(result_RX2);
    [fpr3, tpr3,thre3, auc3, axujing3] = myPlot3DROC(map,result);
    a3_Result(i,1) = auc3;
    a3_Result(i,2) = axujing3;
end
%%
for k = 1:19
    for i = 1:row
        for j = 1:col
            omega(i,j) = 1-exp(-k*O(i,j));
        end
    end
    result_3 = omega.*result;
    output_result = hyperNormalize(result_3);
    a4_Result(k,1) = auc4;
    a4_Result(k,2) = axujing4;
end
a4_maxresult = max(a4_Result);
a4_pos = find(a4_Result==max(a4_Result));
%% 融合
for k = 1:100
    for i = 1:row
        for j = 1:col
            omega(i,j) = 1-exp(-k*initial_det(i,j));
        end
    end
    result_i = gen.*omega;
    a7_Result(k,1) = auc7;
    a7_Result(k,2) = axujing7;
end
auc7 = max(a7_Result);
result=A./(A+exp(-result_i.*Pred_CM));