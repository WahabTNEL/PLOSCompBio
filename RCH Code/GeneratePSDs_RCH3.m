%% Visualize PSDS FOR RCH3 (Generates Supplementary Figure S1)
fs=4;
windowsize=30*fs;
overlap=0.75;
load('Data/RCH3_cell_normalized.mat');
%%
tr=Inf;
[Pxx_all_d,  Power_bands_d, f]=PwelchCell(CellD(1:4), windowsize, overlap, fs, tr);
[Pxx_all_r,  Power_bands_r, f]=PwelchCell(CellR(1:12), windowsize, overlap, fs, tr);
[Pxx_all_el, Power_bands_el,  f]=PwelchCell(CellEL(1:16), windowsize, overlap, fs, tr);
%%
%%
chans=[ 3 8 11 16 28 31 34 38 44 46];

for i=1:10
    
    subplot(2,5,i)
c=chans(i);

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
xlim([-1.8 0.28]);

hold on;
max_val=max([max(mean(Pxx_all_dc,2)) max(mean(Pxx_all_rc,2)) max(mean(Pxx_all_elc,2))]);
min_val=min([min(mean(Pxx_all_dc,2)) min(mean(Pxx_all_rc,2)) min(mean(Pxx_all_elc,2))]);
%plot([log10(1) log10(1)], 10*log10([min_val max_val]), '--','LineWidth', 2, 'Color', 'k');
plot([log10(1/3) log10(1/3)], [-8.5 -2], '--','LineWidth', 2, 'Color', 'k');
%ylim(10*log10([min_val max_val])-[1 3])
ylim([-8.5 -2]);
     
     
     
%legend('Dialogue', 'Rest', 'Electronics');
xlabel('Oscillation Period (s)');
ylabel('dB/Hz');
title(['Channel ' good_channels(c)]);
end
