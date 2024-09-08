clc;
clear all;
close all;
%% HYDICE
load('.\Data\HYDICE_data.mat');
load('.\Data\initial_det.mat');
%% 
I = hyperNormalize(data);
[m, n, z] = size(I);
scale = 0.4;
initial_bw = im2bw(initial_det,0.26);
initial_L = reshape(initial_bw,1,m*n);
initial_S = initial_L(:,1:scale*m*n);
TR = initial_S';
TE = initial_L';
[row,col,bands]=size(data);
I2d = hyperConvert2d(I);
for i = 1 : z
    I2d(i,:) = mat2gray(I2d(i,:));
end
TR2d = hyperConvert2d(TR);
TE2d = hyperConvert2d(TE);
%% 
TE_sample = reshape(data,row*col,bands);
TR_sample = TE_sample(1:scale*m*n,:);
%% 
TE_temp = im2bw(reshape(data,row*col,bands));
TR_temp = im2bw(TE_sample(1:scale*m*n,:));
K = 10;
si = 1;
Train_W = creatLap(TR_sample', K, si);
Train_D = (sum(Train_W, 2)).^(-1/2);
Train_D = diag(Train_D);
L_temp = Train_W * Train_D;
Train_L = Train_D * L_temp;
Train_L = Train_L + eye(size(Train_L));
Test_W = creatLap(TE_sample', K, si);
Test_D = (sum(Test_W, 2)).^(-1/2);
Test_D = diag(Test_D);
L_temp = Test_W * Test_D;
Test_L = Test_D * L_temp;
Test_L = Test_L + eye(size(Test_L));
Train_X = TR_sample;
Test_X = TE_sample;
TrLabel = initial_S';
TeLabel = reshape(initial_bw,row*col,1);
save('.\Input\HYDICE\Train_X.mat','Train_X');
save('.\Input\HYDICE\Test_X.mat','Test_X');
save('.\Input\HYDICE\TrLabel.mat','TrLabel');
save('.\Input\HYDICE\TeLabel.mat','TeLabel');
save('.\Input\HYDICE\Train_L.mat','Train_L');
save('.\Input\HYDICE\Test_L.mat','Test_L');
%% 
load('.\Data\HYDICE_data.mat');
[row,col,bands]=size(data);
data_h2 = data;
%%
for i = 1:32
    data_s=double(data_h2(:,:,i));
    data_s_N = hyperNormalize(data_s);
    Lambda = [50];
    VectorAttributes(1) = 'a';
    NumberOfEAP = 1;
    Lambda = double(sort(nonzeros(Lambda))');
    PC30_int16 = ConvertFromZeroToOneThousand(data_s,false);
    PC30_int16 = int16(PC30_int16);
    AP30 = attribute_profile(PC30_int16,'a', Lambda);
    D(:,:,1) = abs(AP30(:,:,1)-AP30(:,:,2))+abs(AP30(:,:,2)-AP30(:,:,3));
    img = D(:,:,1);

    if length(size(img))>2
        img=rgb2gray(img);
    end
    guild = mat2gray(double(img));
    p =  mat2gray(double(img));
    r = 1;
    eps = 0.3^2;
    O = guidedfilter(guild, p, r, eps);
    O_N = hyperNormalize(O);
    output = im2bw(O_N,0.21);
end
% a1_maxresult = max(a1_Result);
% a1_pos = find(a1_Result==max(a1_Result));
% a2_maxresult = max(a2_Result);
% a2_pos = find(a2_Result==max(a2_Result));
% figure,imagesc(O_N),colormap gray,axis image,axis off,title('cujian');
% figure,imagesc(output),colormap gray,axis image,axis off,title('yuzhi');
%% 选背景和目标,粗检结果百分比
% temp=map;
% temp=output;
% data_1 = reshape(data,row*col,bands);
% temp_1 = reshape(temp,row*col,1);
% data_2 = data_1(temp_1==0,:);
% [m n] = size(data_2);
% data_3 = data_2(1:floor(0.8*m),:);
% allbkg_cj = hyperNormalize(data_2);
% [a,b] = size(allbkg_cj);
% allbkg_cj9 = zeros(floor(0.9*a),b);
% for i = 1:bands
%     x = allbkg_cj(:,i);
%     y = sort(x);
%     alpha = y(floor(0.9*a));
%     z = x(x<=alpha);
%     allbkg_cj9(:,i) = z(1:floor(0.9*a),:);
% end