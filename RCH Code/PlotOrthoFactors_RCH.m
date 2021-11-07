
%% Generate Single Factor performance.

%load('mat_results/run002/gpfa_xDim20.mat');
           
            %%
X=[seqTrain.xsm];
L=estParams.C;
[Xorth, Lorth, DD, TT] = orthogonalize_v2(X, L);
seqTrain=segmentByTrial(seqTrain,Xorth, 'xorth');
trueYtrain=[ones(1,150) 2*ones(1,150) 3*ones(1,150)];

var_exp=round(diag(DD.^2)./sum(diag(DD.^2))*100);
%% Investigate PSDs

windowsize=120;
overlap=0;
fs=4;
tr=Inf;
CellD{1}=[seqTrain(find(trueYtrain==1)).xorth];
CellR{1}=[seqTrain(find(trueYtrain==2)).xorth];
CellEL{1}=[seqTrain(find(trueYtrain==3)).xorth];

%%
[Pxx_all_d,~, f]=PwelchCell(CellD, windowsize, overlap, fs, tr);
[Pxx_all_r, ~,f]=PwelchCell(CellR, windowsize, overlap, fs, tr);
[ Pxx_all_el, ~, f]=PwelchCell(CellEL, windowsize, overlap, fs, tr);
%

%%
for i=1:10               
    
    subplot(2,5,i);

%%

%c=sig_fac(i);
c=i;
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

t=[60 30 15 10 5  3 2 1 0.6];

xticks(log10(1./t));
xticklabels(t);
xlim([-1.8 0.29]);

if c==1
legend('Dialogue', 'Rest', 'Electronics');
end
xlabel('Oscillation Period (s)');
ylabel('dB/Hz');
title(['Factor ' num2str(c) ', Variance Exp= ' num2str(var_exp(i)) '%']);
end