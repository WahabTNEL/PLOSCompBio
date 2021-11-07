
%% GPFA fitting with Cross Validation
%load('D:\Dropbox\Gilja Lab\Second Paper Documents\Code\Context Decoding Code\NYU Code\Data\NYU_epoch_normalized.mat')
%addpath('D:\Dropbox\Gilja Lab\Second Paper Documents\Code\Context Decoding Code');
%addpath('D:\Dropbox\Gilja Lab\Second Paper Documents\Code\Context Decoding Code\RCH Code');
%addpath(genpath('D:\Dropbox\Gilja Lab\Second Paper Documents\Code\Context Decoding Code\gpfa_v0203'));
%% Segment Data according to context
Xd=X(find(Y==0));
Xr=X(find(Y==1));
Xel=X(find(Y==2));
param1=3;
param2=0.5;
param3=1/120;

[bad_epochs_d] =ArtifactRejection(Xd, param1, param2, param3);
[bad_epochs_r] =ArtifactRejection(Xr, param1, param2, param3);
[bad_epochs_el] =ArtifactRejection(Xel, param1, param2, param3);

[Xd_new] =RemoveZeros(Xd);
[Xr_new] =RemoveZeros(Xr);
[Xel_new] =RemoveZeros(Xel);

Xd_new(find(bad_epochs_d==1))=[];
Xr_new(find(bad_epochs_r==1))=[];
Xel_new(find(bad_epochs_el==1))=[];
%%


for i=1:length(Xd_new)
    data_d(i).trialId=i;
    temp=Xd_new{i};
    data_d(i).y=temp;
    data_d(i).T=length(temp);
    
    
end

for i=1:length(Xr_new)
    data_r(i).trialId=i;
    temp=Xr_new{i};
    data_r(i).y=temp;
    data_r(i).T=length(temp);
     
end

for i=1:length(Xel_new)
    data_el(i).trialId=i;
    temp=Xel_new{i};
    data_el(i).y=temp;
    data_el(i).T=length(temp);
end

data=[data_d data_r data_el];
Y=[zeros(1,length(data_d)) ones(1, length(data_r)) 2*ones(1, length(data_el))];
%% Set CV Folds
  K=7;

  
  indD=sort(crossvalind('Kfold', 1:length(find(Y==0)), K));
  indR=sort(crossvalind('Kfold', 1:length(find(Y==1)), K));
  indEL=sort(crossvalind('Kfold', 1:length(find(Y==2)), K));

%%
% ========================================================
% 2) Full cross-validation to find:
%  - optimal state dimensionality for all methods
%  - optimal smoothing kernel width for two-stage methods
% ========================================================

kernSDList=[1 250 500 1000 2000 4000 10000];
% Perform cross-validation for different state dimensionalities.
% Results are saved in mat_results/runXXX/, where XXX is runIdx.
xrange= [2 5 10 15 20];
runIdx=2;
K=7;
parfor i=1:5
  xDim=xrange(i);
 
      neuralTraj_RCH_v2(runIdx, data, indD, indR, indEL,  'method',   'fa', 'xDim', xDim, 'K', K, 'Y', Y , 'emMaxIters', 500, 'binWidth', 250, 'kernSDList', kernSDList);
      neuralTraj_RCH_v2(runIdx, data, indD, indR, indEL,  'method', 'gpfa', 'xDim', xDim, 'K', K, 'Y', Y , 'emMaxIters', 500, 'binWidth', 250, 'kernSDList', kernSDList);
end
fprintf('\n');


%%
plotPredErrorVsDim_v2(1,1);

