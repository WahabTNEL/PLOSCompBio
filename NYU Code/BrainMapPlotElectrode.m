function BrainMapPlotElectrode(brain, x, y, elec, istext, Title)
L=length(elec);            %number of datapoints
data=elec; % create example data set with values ranging from 0 to 3.6
indexValue = 0;     % value for which to set a particular color
topColor = [1 0 0];         % color for maximum data value (red = [1 0 0])
indexColor = [0 0 0];       % color for indexed data value (white = [1 1 1])
bottomcolor = [0 1 0];      % color for minimum data value (blue = [0 0 1])
% Calculate where proportionally indexValue lies between minimum and
% maximum values
%largest = max(max(data));
%smallest = min(min(data));
largest=0.59;
smallest=-0.51;
%largest=max(elec);
%smallest=min(elec);
index = L*abs(indexValue-smallest)/(largest-smallest);
% Create color map ranging from bottom color to index color
% Multipling number of points by 100 adds more resolution
customCMap1 = [linspace(bottomcolor(1),indexColor(1),100*index)',...
            linspace(bottomcolor(2),indexColor(2),100*index)',...
            linspace(bottomcolor(3),indexColor(3),100*index)'];
% Create color map ranging from index color to top color
% Multipling number of points by 100 adds more resolution
customCMap2 = [linspace(indexColor(1),topColor(1),100*(L-index))',...
            linspace(indexColor(2),topColor(2),100*(L-index))',...
            linspace(indexColor(3),topColor(3),100*(L-index))'];
customCMap = [customCMap1;customCMap2];  % Combine colormaps

%%

bc=[1 6 10 23 38 39 40 45 46 51 52 57 58 61 62];
channels=1:64;
channels(bc)=[];


ax1=axes;
imshow(rgb2gray(brain));
YTick=ax1.YTick;
XTick=ax1.XTick;
    if istext==1
            acc=zeros(1,64);
        acc(bc)=nan;
        acc(~isnan(acc))=elec;
        img=flipud(fliplr(reshape(acc, [8 8])'));
    else
    %Using Negative and Positive values
        acc=zeros(1,64);
        acc(bc)=nan;
        elec=round((elec-smallest)/(largest-smallest)*length(customCMap));
        acc(~isnan(acc))=elec;
        img=flipud(fliplr(reshape(acc, [8 8])'));
    end


RGB=ind2rgb(double(img), customCMap);

ax2=axes;

ax2.Visible='off';
ax2.YTick=YTick;
ax2.YDir='reverse';
ax2.Position=ax1.Position;
ax2.XTick=XTick;
linkaxes([ax1, ax2]);

for i=1:64
    [k,j]=ind2sub([8 8], i);
    
    %RGBVec(i,:)=squeeze(RGB(k,j,:))';
    hold on;
    
    %if all(i-bc)==0
     if isnan(img(j,k))==1
        
        scatter(ax2, x(i), y(i), 500, [0 0 0], 'filled');
        
    elseif isnan(img(j,k))==0     
        
        if istext==1
        scatter(ax2, x(i), y(i), 500,  [0 0 1], 'filled');   
        text(x(i),y(i), num2str(img(j,k)), 'Color', [1 1 1],'HorizontalAlignment','center','FontSize',16, 'BackgroundColor', 'none');
        
        else
            scatter(ax2, x(i), y(i), 500,  squeeze(RGB(j,k,:))', 'filled');

        end
        
    end
    

end
if istext==0
    colorbar;
    %
    set([ax1,ax2], 'Position', [0.17 .11 .7 .8]);
    colormap(ax2, customCMap);
    colorbar(ax2, 'Position', [ 0.07 0.2 0.05 0.7], 'Ticks', [0 1],'TickLabels', {num2str(smallest) num2str(largest)});
end
title(ax1,Title);


end






%viscircles(Centroids, ones(64,1)*20);
