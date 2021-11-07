function [data_mean, data_std]=SegmentMeanStd(data, fs, windowsize)
   
    remove=rem(length(data), windowsize*fs);
    
    data_cut=data(:,1:end-remove);
    
    data_segment=reshape(data_cut, [size(data_cut,1)  windowsize*fs length(data_cut)/(windowsize*fs)]);
    
    data_mean=squeeze(mean(data_segment,2));
    data_std=squeeze(std(data_segment, [], 2));
end
