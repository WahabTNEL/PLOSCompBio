%% FOR NYU394 (Subject 1)
% This code fits a GPFA model on the training data and uses the learned
% parameters to classify different behavioral states.
% Use this to geenerate figures 5 and 6.


load('Data/NYU_epoch_normalized.mat')
%% Segment Data according to context
Xd=X(find(Y==0));
Xr=X(find(Y==1));
Xel=X(find(Y==2));
Xtv=X(find(Y==3));


%% Artefact Rejection
param1=3;
param2=0.5;
param3=1/120;

[bad_epochs_d] =ArtifactRejection(Xd, param1, param2, param3);
[bad_epochs_r] =ArtifactRejection(Xr, param1, param2, param3);
[bad_epochs_el] =ArtifactRejection(Xel, param1, param2, param3);
[bad_epochs_tv] =ArtifactRejection(Xtv, param1, param2, param3);

Xd_new=Xd;
Xr_new=Xr;
Xel_new=Xel;
Xtv_new=Xtv;

Xd_new(find(bad_epochs_d==1))=[];
Xr_new(find(bad_epochs_r==1))=[];
Xel_new(find(bad_epochs_el==1))=[];
Xtv_new(find(bad_epochs_tv==1))=[];
T=120;
%% Set up variables for each behavioral context


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


for i=1:length(Xtv_new)
    data_tv(i).trialId=i;
    temp=Xtv_new{i};
    data_tv(i).y=temp;
    data_tv(i).T=length(temp);
end


%% Create CV indices
K=7;
indD=sort(crossvalind('Kfold', 1:length(Xd_new), K));
indR=sort(crossvalind('Kfold', 1:length(Xr_new), K));
indEL=sort(crossvalind('Kfold', 1:length(Xel_new), K));
indTV=sort(crossvalind('Kfold', 1:length(Xtv_new), K));


%% Number of GPFA factors
xDim=30; 

%% Distance of one fold between train and test sets
d=1;

%% Applies cross validation and trains the GPFA model on data, results might slightly differ from manuscript figures
% Due to the random sampling but general trends and conclusions will remain valid.
          
for  k=1:K

            data_d_train=shuffle_trials(data_d(find(indD~=k)), data_type{d});
            data_r_train=shuffle_trials(data_r(find(indR~=k)), data_type{d});
            data_el_train=shuffle_trials(data_el(find(indEL~=k)), data_type{d});
            data_tv_train=shuffle_trials(data_tv(find(indTV~=k)), data_type{d});
            
            minlength=min([length(data_d_train) length(data_r_train)  length(data_el_train)  length(data_tv_train)]);

            data_d_train_b=data_d_train(randperm(length(data_d_train),minlength));
            data_r_train_b=data_r_train(randperm(length(data_r_train),minlength));
            data_el_train_b=data_el_train(randperm(length(data_el_train),minlength));
            data_tv_train_b=data_tv_train(randperm(length(data_tv_train),minlength));
            
            
            data_train_b=[data_d_train_b data_r_train_b data_el_train_b data_tv_train_b];
            trueYtrain=[ones(1, length(data_d_train_b)) 2*ones(1, length(data_r_train_b)) 3*ones(1, length(data_el_train_b)) 4*ones(1, length(data_tv_train_b))];
            %% Generate GPFA model for each context

           [estParams{k}, seq_train]=GPFAtrain(data_train_b,xDim);
           tau=250./sqrt(estParams{k}.gamma);
           [~,idx_fast]=sort(tau,'ascend');
           [~,idx_slow]=sort(tau,'descend');
           tau_fold(:,k)=tau(idx_slow);

           %% find log likelihood for each model
            data_d_test=data_d(find(indD==k));
            data_r_test=data_r(find(indR==k));
            data_el_test=data_el(find(indEL==k));
            data_tv_test=data_tv(find(indTV==k));

            minlength=min([length(data_d_test) length(data_r_test)  length(data_el_test)  length(data_tv_test)]);

            data_d_test_b=data_d_test(randperm(length(data_d_test),minlength));
            data_r_test_b=data_r_test(randperm(length(data_r_test),minlength));
            data_el_test_b=data_el_test(randperm(length(data_el_test),minlength));
            data_tv_test_b=data_tv_test(randperm(length(data_tv_test),minlength));

        
            data_test=shuffle_trials([data_d_test_b data_r_test_b data_el_test_b data_tv_test_b], data_type{dt});
            trueYtest=[ones(1, length(data_d_test_b)) 2*ones(1, length(data_r_test_b)) 3*ones(1, length(data_el_test_b)) 4*ones(1, length(data_tv_test_b))];
            
            [seq_test, LL] = exactInferenceWithLL(data_test, estParams{k});
            %%
            for i=1:xDim
                [~,acc_fast(i,k)]=ClassifyFactors(seq_train, trueYtrain, seq_test, trueYtest, T, idx_fast(1:i), 'quadratic');
                [~,acc_slow(i,k)]=ClassifyFactors(seq_train, trueYtrain, seq_test, trueYtest, T, idx_slow(1:i), 'quadratic');
                [~,acc_fast_linear(i,k)]=ClassifyFactors(seq_train, trueYtrain, seq_test, trueYtest, T, idx_fast(1:i), 'linear');
                [~,acc_slow_linear(i,k)]=ClassifyFactors(seq_train, trueYtrain, seq_test, trueYtest, T, idx_slow(1:i), 'linear');
                [~,acc_slow_sing(i,k)]=ClassifyFactors(seq_train, trueYtrain, seq_test, trueYtest, 120, idx_slow(i), 'quadratic');                

                i
            end

end
        
%% This code generates Figure 5a
options.alpha      = 0.4;
options.line_width = 2;
options.error      = 'sem';


plot_areaerrorbar(acc_slow', options);


options.color_area = [0 0 255]./255;    % Blue theme
options.color_line = [0 0 255]./255;

plot_areaerrorbar(acc_fast', options);
hold on;

options.color_area = [255 0 0]./255;    % Blue theme
options.color_line = [255 0 0]./255;

plot_areaerrorbar(acc_slow', options);
hold on;
options.color_area = [0 255 0]./255;    % Blue theme
options.color_line = [0 255 0]./255;

plot_areaerrorbar(acc_fast_linear', options);
hold on;

options.color_area = [255 131 0]./255;    % Blue theme
options.color_line = [255 131 0]./255;
plot_areaerrorbar(acc_slow_linear', options);
xlim([1 30]);
legend('from fastest (Quadratic)', 'from slowest (Quadratic)', 'from fastest (linear)', 'from slowest (linear)');
xlabel('Number of factors');
ylabel('accuracy');

%% This code generates figure 6a
fig=figure;
left_color = [0 0 1];
right_color = [255 69 0]./255;
set(fig,'defaultAxesColorOrder',[left_color; right_color]);

yyaxis left
bar(mean(acc_slow_sing,2) , 'b');
xlabel('Slowest to Fastest Factors');
ylabel('accuracy');
ylim([0.25 0.35]);
yyaxis right
plot(mean(tau_fold,2), 'Color', right_color, 'LineWidth', 4)
ylabel('tau');
xlim([0 31]);