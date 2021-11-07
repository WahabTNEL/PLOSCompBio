

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
        
    elseif img(j,k)>0     
        
        if istext==1
        scatter(ax2, x(i), y(i), 500,  [0 0 1], 'filled');   
        text(x(i),y(i), num2str(img(j,k)), 'Color', [1 1 1],'HorizontalAlignment','center','FontSize',16, 'BackgroundColor', 'none');
        
        else
            scatter(ax2, x(i), y(i), 500,  squeeze(RGB(j,k,:))', 'filled');

        end
        
    end
    

end



