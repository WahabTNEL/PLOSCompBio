%% This code visualizes the factor weights generated from the GPFA model
    
   % load('~/mat_results/run001/gpfa_xDim30.mat');
    tau=250./sqrt(estParams.gamma);
    [~,idx_fast]=sort(tau,'ascend');
    [~,idx_slow]=sort(tau,'descend');
    figure;
    imagesc(estParams.C(:,idx_slow));
    [customCMap]=CreateColorMap(estParams.C);
    colormap(customCMap);    colorbar;
    load('ElectrodeCentroidLocs.mat')
    xlabel('Factors (Slowest to Fastest)');

    %%
    %load('~/mat_results/run002/gpfa_xDim20.mat');
    figure;
    tau=250./sqrt(estParams.gamma);
    [~,idx_fast]=sort(tau,'ascend');
    [~,idx_slow]=sort(tau,'descend');
    imagesc(estParams.C(:,idx_slow));
    %imagesc(Lorth(:,idx_slow));
    [customCMap]=CreateColorMap(estParams.C);
    colormap(customCMap);    
    colorbar;
    yticks([1:26]);
    yticklabels(good_channels);
    xlabel('Factors (Slowest to Fastest)');
    ylabel('Channel');