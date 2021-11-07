%% This code plots the fastest and slowest factors for NYU
%load('/home/wahab/Desktop/Main Code Base/mat_results/run001/gpfa_xDim30.mat');

%%
   tau=250./sqrt(estParams.gamma);
   [~,idx_fast]=sort(tau,'ascend');
    [~,idx_slow]=sort(tau,'descend');
    figure;
    imagesc(estParams.C(:,idx_slow));
    [customCMap]=CreateColorMap(estParams.C);
    colormap(customCMap);    colorbar;
    load('ElectrodeCentroidLocs.mat')
    xlabel('Factors (Slowest to Fastest)');
    ylabel('Channel Number');
    brain=imread('NY394_T1_GS_lateral_rh.png');
    %brain=imread('P1rightw.jpg');
      %%
e=101;
    elec=estParams.C(:,idx_slow(1));
    figure;
    BrainMapPlotElectrode(brain, x, y,  elec, 0, '')
    figure;
    t=0.25:30/120:30;
    plot(t,seqTrain(e).xsm(idx_slow(1),:), 'LineWidth', 2)
    hold on;
    plot(t,seqTrain(e).y(39,:), 'r', 'LineWidth', 1.5);
    xlabel('Time (s)');
    ylabel('Amplitude a.u');
    title('Slowest latent dimension');
    legend('Latent Factor', 'Single Channel');
set(gca, 'FontSize', 18);
%%

e=1;
    elec=estParams.C(:,idx_slow(end));
    figure;
    BrainMapPlotElectrode(brain, x, y,  elec, 0, '')
    figure;
    t=0.25:30/120:30;
    plot(t,seqTrain(e).xsm(idx_slow(end),:)*-1, 'LineWidth', 2)
    hold on;
    plot(t,seqTrain(e).y(3,:), 'g', 'LineWidth', 1.5);
    xlabel('Time (s)');
    ylabel('Amplitude a.u');
    title('Fastest latent dimension');
    legend('Latent Factor', 'Single Channel');
set(gca, 'FontSize', 18);



    