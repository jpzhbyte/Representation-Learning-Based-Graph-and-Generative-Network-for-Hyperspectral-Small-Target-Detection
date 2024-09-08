%     demo for CSCR detection algorithm
%--------------Brief description-------------------------------------------
%
% 
% This demo implements the CSCR hyperspectral image target detection [1]
% 
%
% More details in:
%
% [1] W. Li, Q. Du, and B. Zhang, Combined Sparse and
% Collaborative Representation for Hyperspectral Target Detection.
% Pattern Recognition, vol. 48, no. 12, pp. 3904-3916, December 2015. 
% 
%
% contact: liwei089@ieee.org (liw@mail.buct.edu.cn) (Wei Li)


clear all;clc;close all;

d = d';
[a b c] = size(data);
map_1=zeros(a,b);

Data_Ori=data(:,:,:);
d=d';
TargetTrain=d(:,:)';

class_num = 6;
for i=1:1:a
    for j=1:1:b
        if(map(i,j)==class_num)
        map_1(i,j)=1;
        end
    end
end

tic
lambda = 1e-2; beta = 1e-2;
output = CR_SUNSAL_Detection_Local(Data_Ori, TargetTrain, 15, 5, lambda, beta); 
output3 = reshape(output, a, b, 2);
tmp1 = output3(:,:,1); tmp2 = output3(:,:,2);
output3(:,:,1) = (tmp1./max(tmp1(:)));
output3(:,:,2) = (tmp2./max(tmp2(:)));

toc
m=output3(:,:,1);
for i=1:1:a
    for j=1:1:b
        if(isnan(output3(i,j,1))==1)
            output3(i,j,1)=0;
        end
    end
end
r1 = reshape(output3(:,:,1), 1, a*b);

r=reshape(r1,a,b);
sigma1 = 7;
gausFilter = fspecial('gaussian',[5 5],sigma1);
r=imfilter(r,gausFilter,'replicate');
cscr_pa9=r;
save('pavia_cscr9','cscr_pa9');
figure(2)
imagesc(mat2gray(r))
colormap gray,axis off,axis image;
figure(3)
imagesc(mat2gray(map))
colormap gray,axis off,axis image;

N=a*b;
targets = reshape(map, 1, N); 
outputs = reshape(mat2gray(r), 1, N);
[FPR,TPR,thre,auc] = myPlot3DROC(targets, outputs);
cscr = reshape(outputs,a,b);
auc =-trapz(FPR,TPR);
fpr = -trapz(FPR,thre);
figure,plot(FPR,TPR);
xlabel('false alarm rate');
ylabel('probability of detection');
title('ROC curves of detection algorithm');
text(0.2,0.7,sprintf('AUCROC_ =%f',auc));