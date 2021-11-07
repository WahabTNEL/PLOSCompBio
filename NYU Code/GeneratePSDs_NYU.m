%% This code is to visualize the power spectral densities of the band limited amplitude for NYU394 (Subject 1) 
% Plots Figure 2a

%We specify the new sampling rate = 4 Hz, window size for the pwelch method
%of PSD estimation and the amount of overlap
fs=4;
windowsize=30*fs;
overlap=0.75;

%Load data (high gamma amplitude estimate that is averaged across 250 ms
%bins)
load('Data/NYU_cell_normalized.mat');
cd ..
addpath('NYU Code');
%%
%Specify outlier rejection parameter tr
tr=5;

%Generate PSDs from the binned high gamma amplitude estimates of the 4
%behavioral states. Each continuous behavioral state is stored in a cell.
[Pxx_all_d,  Power_bands_d, f]=PwelchCell(CellD(1:24), windowsize, overlap, fs, tr);
[Pxx_all_r,  Power_bands_r, f]=PwelchCell(CellR(1:30), windowsize, overlap, fs, tr);
[Pxx_all_el,  Power_bands_el, f]=PwelchCell(CellEL(1:8), windowsize, overlap, fs, tr);
[Pxx_all_tv,  Power_bands_tv, f]=PwelchCell(CellTV(1:12), windowsize, overlap, fs,tr);


%% Plot the PSDs for an examplar channel in log-log scale, the x-axis has been changed to oscillation period (s) rather than frequency (Hz)
c=33;

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
plot([log10(1/3) log10(1/3)], [-12 6], 'k', 'LineStyle', '--', 'LineWidth', 2)
legend('Dialogue', 'Rest', 'Electronics', 'TV');
xlabel('Oscillation Period (s)');
ylabel('dB/Hz');
load('Data/NYU_epoch_normalized', 'channel_idx');
title(['Channel ' num2str(channel_idx(c))]);