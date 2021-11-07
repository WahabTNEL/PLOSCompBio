%%For NYU394 Subject 1
%This code takes the normalized data, estimates the PSD, and using the
%power at slow and fast frequency bands as features in a classification
%scheme for examplar channels for visualization purposes.

load('Data/NYU_cell_normalized.mat');
cd ..
addpath('NYU Code');
%% Normalize Data to remove mean and variance effects
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

TV_low=[];
TV_high=[];

C=size(CellD{1},1);

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

for i=1:length(CellTV)
    
   temp=CellTV{i};
    TV_low1=abs(hilbert(filter(b_low, a_low, temp, [], 2)'))';
    TV_high1=abs(hilbert(filter(b_high, a_high, temp, [], 2)'))';
    [TV_low_cell, ~]=SegmentMeanStd(TV_low1, fs, windowsize);
    [TV_high_cell, ~]=SegmentMeanStd(TV_high1, fs, windowsize);
    TV_low=[TV_low TV_low_cell];
    TV_high=[TV_high TV_high_cell];
    
end
%%
K=7;
d=1;
numofFactors=1;
indD=sort(crossvalind('Kfold', squeeze(D_low(1,:)), K));
indR=sort(crossvalind('Kfold', squeeze(R_low(1,:)), K));
indEL=sort(crossvalind('Kfold', squeeze(EL_low(1,:)), K));
indTV=sort(crossvalind('Kfold', squeeze(TV_low(1,:)), K));
fprintf('Classifying mean...  \n');

%%
fprintf('Classifying both...  \n');

for i=1:K
  
    for c=1:C
        
            DTrain1=D_high(c,  find(indD<i-d | indD>i+d));
            RTrain1=R_high(c, find(indR<i-d | indR>i+d));
            ELTrain1=EL_high(c,  find(indEL<i-d | indEL>i+d));
            TVTrain1=TV_high(c,  find(indTV<i-d | indTV>i+d));
            %
            minlength=min([size(DTrain1,2) size(RTrain1,2) size(ELTrain1,2) size(TVTrain1,2)]);

            DTrain_h=DTrain1(randperm(size(DTrain1,2),minlength));
            RTrain_h=RTrain1(randperm(size(RTrain1,2),minlength));
            ELTrain_h=ELTrain1(randperm(size(ELTrain1,2),minlength));
            TVTrain_h=TVTrain1(randperm(size(TVTrain1,2),minlength));



            DTest1=D_high(c,find(indD==i));
            RTest1=R_high(c,find(indR==i));
            ELTest1=EL_high(c,find(indEL==i));
            TVTest1=TV_high(c,find(indTV==i));


            minlength=min([size(DTest1,2) size(RTest1,2) size(ELTest1,2) size(TVTest1,2)]);

            DTest_h=DTest1(randperm(size(DTest1,2),minlength));
            RTest_h=RTest1(randperm(size(RTest1,2),minlength));
            ELTest_h=ELTest1(randperm(size(ELTest1,2),minlength));
            TVTest_h=TVTest1(randperm(size(TVTest1,2),minlength));
            
            DTrain1=D_low(c,  find(indD<i-d | indD>i+d));
            RTrain1=R_low(c, find(indR<i-d | indR>i+d));
            ELTrain1=EL_low(c,  find(indEL<i-d | indEL>i+d));
            TVTrain1=TV_low(c,  find(indTV<i-d | indTV>i+d));
            %
            minlength=min([size(DTrain1,2) size(RTrain1,2) size(ELTrain1,2) size(TVTrain1,2)]);

            DTrain_l=DTrain1(randperm(size(DTrain1,2),minlength));
            RTrain_l=RTrain1(randperm(size(RTrain1,2),minlength));
            ELTrain_l=ELTrain1(randperm(size(ELTrain1,2),minlength));
            TVTrain_l=TVTrain1(randperm(size(TVTrain1,2),minlength));



            DTest1=D_low(c,find(indD==i));
            RTest1=R_low(c,find(indR==i));
            ELTest1=EL_low(c,find(indEL==i));
            TVTest1=TV_low(c,find(indTV==i));


            minlength=min([size(DTest1,2) size(RTest1,2) size(ELTest1,2) size(TVTest1,2)]);

            DTest_l=DTest1(randperm(size(DTest1,2),minlength));
            RTest_l=RTest1(randperm(size(RTest1,2),minlength));
            ELTest_l=ELTest1(randperm(size(ELTest1,2),minlength));
            TVTest_l=TVTest1(randperm(size(TVTest1,2),minlength));
            
            DTrain=[DTrain_h; DTrain_l];
            RTrain=[RTrain_h; RTrain_l];
            ELTrain=[ELTrain_h; ELTrain_l];
            TVTrain=[TVTrain_h; TVTrain_l];
            
            DTest=[DTest_h; DTest_l];
            RTest=[RTest_h; RTest_l];
            ELTest=[ELTest_h; ELTest_l];
            TVTest=[TVTest_h; TVTest_l];
            
            [~, accuracyM4Context_both(i,c)]=ClassifyFAMultiClass(DTrain, RTrain, ELTrain, TVTrain, DTest, RTest, ELTest, TVTest, numofFactors, 'No');
            [~, accuracyM4Context_high(i,c)]=ClassifyFAMultiClass(DTrain_h, RTrain_h, ELTrain_h, TVTrain_h, DTest_h, RTest_h, ELTest_h, TVTest_h, numofFactors, 'No');
            [~, accuracyM4Context_low(i,c)]=ClassifyFAMultiClass(DTrain_l, RTrain_l, ELTrain_l, TVTrain_l, DTest_l, RTest_l, ELTest_l, TVTest_l, numofFactors, 'No');
    end
    
i   
end

%% Visualize PSDS
fs=4;
windowsize=30*fs;
%windowsize=257;
overlap=0.75;
load('Data/NYU_cell_normalized.mat');
%%
tr=5;
[Pxx_all_d,  Power_bands_d, f]=PwelchCell(CellD(1:24), windowsize, overlap, fs, tr);
[Pxx_all_r,  Power_bands_r, f]=PwelchCell(CellR(1:30), windowsize, overlap, fs, tr);
[Pxx_all_el,  Power_bands_el, f]=PwelchCell(CellEL(1:8), windowsize, overlap, fs, tr);
[Pxx_all_tv,  Power_bands_tv, f]=PwelchCell(CellTV(1:12), windowsize, overlap, fs,tr);



%%
[~, idx]=sort(mean(accuracyM4Context_both), 'descend');


for i=1:10
    
        h=subplot(2,5,i);
        c=idx(i);

        options.alpha      = 0.4;
        options.line_width = 2;
        options.error      = 'sem';

        Pxx_all_dc=Pxx_all_d{c};

        options.x_axis=log10(f(2:end));

        options.color_area = [0 0 255]./255;    % Blue theme
        options.color_line = [0 0 255]./255;

        plot_areaerrorbar(10*log10(Pxx_all_dc(2:end,:)'), options)



        options.color_area = [255 0 0 ]./255;    % Orange theme
        options.color_line = [255 0 0]./255;
        hold on;      
        Pxx_all_rc=Pxx_all_r{c};

        options.x_axis=log10(f(2:end));
        plot_areaerrorbar(10*log10(Pxx_all_rc(2:end,:)'), options)       




        options.color_area = [0 255 0 ]./255;    % Orange theme
        options.color_line = [0 255 0]./255;
        hold on;      
        Pxx_all_elc=Pxx_all_el{c};

        options.x_axis=log10(f(2:end));
        plot_areaerrorbar(10*log10(Pxx_all_elc(2:end,:)'), options)

        options.color_area = [255 131 0]./255;    % Orange theme
        options.color_line = [255 131 0]./255;
        hold on;      
        Pxx_all_tvc=Pxx_all_tv{c};

        options.x_axis=log10(f(2:end));
        plot_areaerrorbar(10*log10(Pxx_all_tvc(2:end,:)'), options)      

        t=[60 30 15 10 5  3 2 1 0.6];



        xticks(log10(1./t));
        xticklabels(t);
        xlim([-1.8 0.28]);
        hold on;
        plot([log10(1/3) log10(1/3)], [-12 8], 'k', 'LineStyle', '--', 'LineWidth', 2)
        if i==1
        legend('Dialogue', 'Rest', 'Electronics', 'TV');
        end
        xlabel('Oscillation Period (s)');
        ylabel('dB/Hz');
        load('Data/NYU_epoch_normalized', 'channel_idx');
        title(['Channel ' num2str(channel_idx(idx(i)))]);

        dim = h.Position;
        str = {['Slow = ' num2str(round(mean(accuracyM4Context_low(:,c))*100)) ...
                 newline ' Fast = ' num2str(round(mean(accuracyM4Context_high(:,c))*100)) ...
                 newline ' Both = ' num2str(round(mean(accuracyM4Context_both(:,c))*100))]};
        annotation('textbox',dim,'String',str,'FitBoxToText','on');
        ylim([-11 8]);
end

