function segmented_signal=segment(signal, windowsize, overlap)

    
        N=length(signal);
    
        overlapsize=windowsize-windowsize*overlap;

        
        N_seg=ceil(N/(windowsize-overlap*windowsize));

   
    for i=1:N_seg
        start=(i-1)*overlapsize+1;
        endpoint=start+windowsize-1;
        
        if endpoint>N
            continue
        end
        
        segmented_signal(:,i)=signal(start:endpoint);
    end
    
    
end
