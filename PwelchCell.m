%% This function estimates the PSD from the binned high gamma amplitude estimate.
% This includes some outlier rejection, estimation of the
% power in 3 frequency bands by integrating the PSD

%note that Pxx_dist is the periodogram estimate for each window, to
%estimate the PSD take the average of Pxx_dist


function [Pxx_dist, Power_bands_c, f]=PwelchCell(Data, windowsize, overlap, fs, tr)
    
    ch=size(Data{1},1);

                             
    for c=1:ch
      Pxx_all=[];
      k=1;
         for i=1:length(Data)
           
       
            
            temp=Data{i}(c,:);
           
            temp_seg=segment(temp,windowsize, overlap);
            
            remove=find(isoutlier(mean(temp_seg), 'median', 'threshold', tr));

            temp_seg(:,find(remove))=[];

            
            if ~isempty(temp_seg)==1
                [Pxx1, f]=periodogram(temp_seg,hanning(windowsize), [] , fs);
                ind_S=find(f<=0.2);
                ind_M=find(f>0.2 & f<=1);
                ind_F=find(f>1);
                for j=1:size(Pxx1,2)

                Power_bands(k,1)=trapz(f(ind_S), Pxx1(ind_S,j));
                Power_bands(k,2)=trapz(f(ind_M), Pxx1(ind_M,j));
                Power_bands(k,3)=trapz(f(ind_F), Pxx1(ind_F,j));
        
                 k=k+1;
                end
                 Pxx_all=[Pxx_all Pxx1];
            end
            
            
         end
         

            Pxx_dist{c}=Pxx_all;
            Power_bands_c{c}=Power_bands;
    end
    

    
    
    
    
    
    
end

