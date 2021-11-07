%% FOR RCH1 (Subject 2)
% This code takes the unnormalized binned high gamma amplitude and classifies each 30 second epoch using the mean and variance as features
% The covariance matrix is also estimated in each 30 s epoch from the normalized data and used as a
% feature in another classifier as well as the slow (slower than 1/3 Hz)
% and fast (faster than 1/3 Hz) spectral features of the high gamma
% amplitude estimate.

% Note: due to the random sampling inherent in this method, the final plot
% might not match the figure in the manuscript exactly, but the general
% trends should be equivalent
%Generates Figure 2b
clear all;
addpath('RCH Code');
addpath(genpath('Covariance Toolbox'));
load('RCH/NYU_cell_unnormalized.mat');

%%
fs=4;
trial_length=30;

D_mean=[];
D_std=[];

R_mean=[];
R_std=[];

EL_mean=[];
EL_std=[];

for i=1:length(CellD)
    temp=CellD{i};
    [D_mean_cell, D_std_cell]=SegmentMeanStd(temp, fs, trial_length);
    D_mean=[D_mean D_mean_cell];
    D_std=[D_std D_std_cell];
end
for i=1:length(CellR)
  temp=CellR{i};
  [R_mean_cell, R_std_cell]=SegmentMeanStd(temp, fs, trial_length);
  R_mean=[R_mean R_mean_cell];
  R_std=[R_std R_std_cell];
end
for i=1:length(CellEL)
  temp=CellEL{i};
  [EL_mean_cell, EL_std_cell]=SegmentMeanStd(temp, fs, trial_length);
  EL_mean=[EL_mean EL_mean_cell];
  EL_std=[EL_std EL_std_cell];
end


%% Create CV fold indices
K=7;
indD=sort(crossvalind('Kfold', squeeze(D_mean(1,:)), K));
indR=sort(crossvalind('Kfold', squeeze(R_mean(1,:)), K));
indEL=sort(crossvalind('Kfold', squeeze(EL_mean(1,:)), K));
fprintf('Classifying mean...  \n');
%%
numofFactors=1;
d=1;
for i=1:K
 
   
DTrain1=D_mean(:,  find(indD<i-d | indD>i+d));
RTrain1=R_mean(:,  find(indR<i-d | indR>i+d));
ELTrain1=EL_mean(:,  find(indEL<i-d | indEL>i+d));

minlength=min([size(DTrain1,2) size(RTrain1,2) size(ELTrain1,2)]);

DTrain=DTrain1(:, randperm(size(DTrain1,2),minlength));
RTrain=RTrain1(:, randperm(size(RTrain1,2),minlength));
ELTrain=ELTrain1(:, randperm(size(ELTrain1,2),minlength));
%%
%{  
DTrain=Dhil(:,  find(indD~=i));
RTrain=Rhil(:,  find(indR~=i));
ELTrain=ELhil(:,  find(indEL~=i));
%}
DTest1=D_mean(:, find(indD==i));
RTest1=R_mean(:, find(indR==i));
ELTest1=EL_mean(:, find(indEL==i));


minlength=min([size(DTest1,2) size(RTest1,2) size(ELTest1,2)]);

DTest=DTest1(:, randperm(size(DTest1,2),minlength));
RTest=RTest1(:, randperm(size(RTest1,2),minlength));
ELTest=ELTest1(:, randperm(size(ELTest1,2),minlength));


[CM3Context_mean(:,:,i), accuracyM3Context_mean(i)]=ClassifyFAMThreeClass(DTrain, RTrain, ELTrain,  DTest,  RTest,  ELTest, numofFactors, 'No');
%[CMDEL(:,:,i), accuracyDEL(i)]=BinaryClassifyFA(DTrain, ELTrain,  DTest,  ELTest, numofFactors, 'Yes');
i   
end

%% Test for shifts in variance only

fprintf('Classifying variance...  \n');



for i=1:K
 
   
DTrain1=D_std(:,  find(indD<i-d | indD>i+d));
RTrain1=R_std(:,  find(indR<i-d | indR>i+d));
ELTrain1=EL_std(:,  find(indEL<i-d | indEL>i+d));

minlength=min([size(DTrain1,2) size(RTrain1,2) size(ELTrain1,2)]);

DTrain=DTrain1(:, randperm(size(DTrain1,2),minlength));
RTrain=RTrain1(:, randperm(size(RTrain1,2),minlength));
ELTrain=ELTrain1(:, randperm(size(ELTrain1,2),minlength));
%%
%{  
DTrain=Dhil(:,  find(indD~=i));
RTrain=Rhil(:,  find(indR~=i));
ELTrain=ELhil(:,  find(indEL~=i));
%}
DTest1=D_std(:, find(indD==i));
RTest1=R_std(:, find(indR==i));
ELTest1=EL_std(:, find(indEL==i));


minlength=min([size(DTest1,2) size(RTest1,2) size(ELTest1,2)]);

DTest=DTest1(:, randperm(size(DTest1,2),minlength));
RTest=RTest1(:, randperm(size(RTest1,2),minlength));
ELTest=ELTest1(:, randperm(size(ELTest1,2),minlength));


[CM3Context_std(:,:,i), accuracyM3Context_std(i)]=ClassifyFAMThreeClass(DTrain, RTrain, ELTrain,  DTest,  RTest,  ELTest, numofFactors, 'No');
%[CMDEL(:,:,i), accuracyDEL(i)]=BinaryClassifyFA(DTrain, ELTrain,  DTest,  ELTest, numofFactors, 'Yes');
i   
end

%% Generate Covariance Matrices using normalized high gamma band amplitude estimates and classify

load('Data/RCH1_cell_normalized.mat');

fs=4;
windowsize=30;
c=1;
for i=1:length(CellD)
    temp=CellD{i};
    remove=rem(length(temp), windowsize*fs);
    D_cut=temp(:,1:end-remove);    
    D_segment=reshape(D_cut, [size(D_cut,1)  windowsize*fs length(D_cut)/(windowsize*fs)]);
    for k=1:size(D_segment,3)
        D_sample=squeeze(D_segment(:,:,k));
        Cov_D(:,:,c)=cov(D_sample');
        c=c+1;
    end
end

c=1;
for i=1:length(CellR)
    temp=CellR{i};
    remove=rem(length(temp), windowsize*fs);
    R_cut=temp(:,1:end-remove);    
    R_segment=reshape(R_cut, [size(R_cut,1)  windowsize*fs length(R_cut)/(windowsize*fs)]);
    for k=1:size(R_segment,3)
        R_sample=squeeze(R_segment(:,:,k));
        Cov_R(:,:,c)=cov(R_sample');
        c=c+1;
    end
end

c=1;
for i=1:length(CellEL)
    temp=CellEL{i};
    remove=rem(length(temp), windowsize*fs);
    EL_cut=temp(:,1:end-remove);    
    EL_segment=reshape(EL_cut, [size(EL_cut,1)  windowsize*fs length(EL_cut)/(windowsize*fs)]);
    for k=1:size(EL_segment,3)
        EL_sample=squeeze(EL_segment(:,:,k));
        Cov_EL(:,:,c)=cov(EL_sample');
        c=c+1;
    end
end



%%
addpath(genpath('/home/wahab/Desktop/Main Code Base/Covariance Toolbox'));
d=1;
metric_mean = {'euclid','logeuclid','riemann','ld'};
metric_dist = {'euclid','logeuclid','riemann','ld','kullback'};
acc = zeros(length(metric_mean),length(metric_dist));
for k=1:K
    
DTrain1=Cov_D(:, :, find(indD<k-d | indD>k+d));
RTrain1=Cov_R(:, :, find(indR<k-d | indR>k+d));
ELTrain1=Cov_EL(:, :, find(indEL<k-d | indEL>k+d));

minlength=min([size(DTrain1,3) size(RTrain1,3) size(ELTrain1,3)]);

DTrain=DTrain1(:, :, randperm(size(DTrain1,3),minlength));
RTrain=RTrain1(:, :, randperm(size(RTrain1,3),minlength));
ELTrain=ELTrain1(:, :, randperm(size(ELTrain1,3),minlength));



    COVtrain=cat(3,DTrain, RTrain, ELTrain);
    Ytrain=[ones(1,size(DTrain,3)) 2*ones(1, size(RTrain,3)) 3*ones(1,size(ELTrain,3))];
  
    
DTest1=Cov_D(:, :, find(indD==k));
RTest1=Cov_R(:, :, find(indR==k));
ELTest1=Cov_EL(:,:,  find(indEL==k));


minlength=min([size(DTest1,3) size(RTest1,3) size(ELTest1,3)]);

DTest=DTest1(:, :, randperm(size(DTest1,3),minlength));
RTest=RTest1(:, :, randperm(size(RTest1,3),minlength));
ELTest=ELTest1(:, :, randperm(size(ELTest1,3),minlength));

% Data formating
   COVtest = cat(3, DTest, RTest, ELTest);
   trueYtest  = [ones(1,size(DTest,3)) 2*ones(1, size(RTest,3)) 3*ones(1,size(ELTest,3))];


%% MDM classification - Multiclass


for i=1:length(metric_mean)
   for j=1:length(metric_dist)

        Ytest = mdm(COVtest,COVtrain,Ytrain,metric_mean{i},metric_dist{j});
        acc(i,j,k) = 100*mean(Ytest==trueYtest);
        [C, order] = confusionmat(trueYtest, Ytest);
        C(1,:) = C(1,:)./sum(C(1,:));
        C(2,:) = C(2,:)./sum(C(2,:));
        C(3,:) = C(3,:)./sum(C(3,:));
        C_cov(:,:,i,j,k)=C;
   end
end

disp('------------------------------------------------------------------');
disp('Accuracy (%) - Rows : distance metric, Colums : mean metric');
disp('------------------------------------------------------------------');
displaytable(squeeze(acc(:,:,k))',metric_mean,10,{'.1f'},metric_dist)
disp('------------------------------------------------------------------');
end


%% Filter the high gamma band amplitude estimate using a lowpass and highpass filter with cutoff 1/Hz,
% Use resulting signal as a feature in an SVM classifier for each 30 second epoch.
D_all=[];
R_all=[];
EL_all=[];

w1=250/1000;
fs=1/w1;
fc=0.333;
N=6;
fs=4;
windowsize=30;

[b_low,a_low] = ellip(N,2,40, fc/(fs/2));
[b_high,a_high] = ellip(N,2,40, fc/(fs/2), 'high');

D_low=[];
D_high=[];

R_low=[];
R_high=[];

EL_low=[];
EL_high=[];

for i=1:length(CellD)
    
    temp=CellD{i};    
    D_low1=abs(hilbert(filter(b_low, a_low, temp, [], 2)'))';
    D_high1=abs(hilbert(filter(b_high, a_high, temp, [], 2)'))';
    [D_low_cell, ~]=SegmentMeanStd(D_low1, fs, windowsize);
    [D_high_cell, ~]=SegmentMeanStd(D_high1, fs, windowsize);
    D_low=[D_low D_low_cell];
    D_high=[D_high D_high_cell];
    
end
for i=1:length(CellR)
    
    temp=CellR{i};
    R_low1=abs(hilbert(filter(b_low, a_low, temp, [], 2)'))';
    R_high1=abs(hilbert(filter(b_high, a_high, temp, [], 2)'))';
    [R_low_cell, ~]=SegmentMeanStd(R_low1, fs, windowsize);
    [R_high_cell, ~]=SegmentMeanStd(R_high1, fs, windowsize);
    R_low=[R_low R_low_cell];
    R_high=[R_high R_high_cell];
    
end
for i=1:length(CellEL)
    
   temp=CellEL{i};
    EL_low1=abs(hilbert(filter(b_low, a_low, temp, [], 2)'))';
    EL_high1=abs(hilbert(filter(b_high, a_high, temp, [], 2)'))';
    [EL_low_cell, ~]=SegmentMeanStd(EL_low1, fs, windowsize);
    [EL_high_cell, ~]=SegmentMeanStd(EL_high1, fs, windowsize);
    EL_low=[EL_low EL_low_cell];
    EL_high=[EL_high EL_high_cell];
    
end




%%
fprintf('Classifying slow...  \n');
numofFactors=10;
for i=1:K
 
   
DTrain1=D_low(:,  find(indD<i-d | indD>i+d));
RTrain1=R_low(:, find(indR<i-d | indR>i+d));
ELTrain1=EL_low(:,  find(indEL<i-d | indEL>i+d));

minlength=min([size(DTrain1,2) size(RTrain1,2) size(ELTrain1,2)]);

DTrain=DTrain1(:, randperm(size(DTrain1,2),minlength));
RTrain=RTrain1(:, randperm(size(RTrain1,2),minlength));
ELTrain=ELTrain1(:, randperm(size(ELTrain1,2),minlength));

DTest1=D_low(:, find(indD==i));
RTest1=R_low(:, find(indR==i));
ELTest1=EL_low(:, find(indEL==i));


minlength=min([size(DTest1,2) size(RTest1,2) size(ELTest1,2)]);

DTest=DTest1(:, randperm(size(DTest1,2),minlength));
RTest=RTest1(:, randperm(size(RTest1,2),minlength));
ELTest=ELTest1(:, randperm(size(ELTest1,2),minlength));


[CM3Context_low(:,:,i), accuracyM3Context_low(i)]=ClassifyFAMThreeClass(DTrain, RTrain, ELTrain,  DTest,  RTest,  ELTest, numofFactors, 'No');
%[CMDEL(:,:,i), accuracyDEL(i)]=BinaryClassifyFA(DTrain, ELTrain,  DTest,  ELTest, numofFactors, 'Yes');
i   
end

%%
fprintf('Classifying fast...  \n');
%%
for i=1:K
 
   
DTrain1=D_high(:,  find(indD~=i));
RTrain1=R_high(:,  find(indR~=i));
ELTrain1=EL_high(:,  find(indEL~=i));

minlength=min([size(DTrain1,2) size(RTrain1,2) size(ELTrain1,2)]);

DTrain=DTrain1(:, randperm(size(DTrain1,2),minlength));
RTrain=RTrain1(:, randperm(size(RTrain1,2),minlength));
ELTrain=ELTrain1(:, randperm(size(ELTrain1,2),minlength));

%{  
DTrain=Dhil(:,  find(indD~=i));
RTrain=Rhil(:,  find(indR~=i));
ELTrain=ELhil(:,  find(indEL~=i));
%}
DTest1=D_high(:, find(indD==i));
RTest1=R_high(:, find(indR==i));
ELTest1=EL_high(:, find(indEL==i));


minlength=min([size(DTest1,2) size(RTest1,2) size(ELTest1,2)]);

DTest=DTest1(:, randperm(size(DTest1,2),minlength));
RTest=RTest1(:, randperm(size(RTest1,2),minlength));
ELTest=ELTest1(:, randperm(size(ELTest1,2),minlength));


[CM3Context_high(:,:,i), accuracyM3Context_high(i)]=ClassifyFAMThreeClass(DTrain, RTrain, ELTrain,  DTest,  RTest,  ELTest, numofFactors, 'No');
%[CMDEL(:,:,i), accuracyDEL(i)]=BinaryClassifyFA(DTrain, ELTrain,  DTest,  ELTest, numofFactors, 'Yes');
i   
end

%% Generate Bar plot

mean_mean=mean(accuracyM3Context_mean,2);
std_mean=std(accuracyM3Context_std,[],2)./sqrt(K);

mean_std=mean(accuracyM3Context_std,2);
std_std=std(accuracyM3Context_std,[],2)./sqrt(K);


accuracyM3Context_cov=squeeze(acc(1,1,:))'/100;

mean_cov=mean(accuracyM3Context_cov,2);
std_cov=std(accuracyM3Context_cov,[],2)./sqrt(K);


mean_low=mean(accuracyM3Context_low,2);
std_low=std(accuracyM3Context_low,[],2)./sqrt(K);

mean_high=mean(accuracyM3Context_high,2);
std_high=std(accuracyM3Context_high,[],2)./sqrt(K);


mean_all=[mean_mean mean_std mean_cov mean_low mean_high];
error_low=[std_mean std_std std_cov std_low std_high];
error_high=[std_mean std_std std_cov std_low std_high]*-1;

x=1:5;

bar(x,mean_all, 'b')                

hold on

er = errorbar(x,mean_all, error_low, error_high, 'LineWidth', 2);    
er.Color = [0 0 0];                            
er.LineStyle = 'none';

ylabel('Accuracy');
%xlabel('Feature');
%title('RCH1');
ylim ([ 0 1])
xticklabels({'Mean', 'Variance', 'Covariance', 'Slow', 'Fast'});
xtickangle(70);
set(gca, 'FontSize', 18);
hold on;
plot([0 6], [0.333 0.333], '--', 'LineWidth', 4)