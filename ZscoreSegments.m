function [z_scored_signal]=ZscoreSegments(signal, windowsize)
    
    % This code Z-scores each segment separately
       z_scored_signal=[];
    
       N=length(signal);
    
   
        
        N_seg=ceil(N/windowsize);

   
    for i=1:N_seg
        start=(i-1)*windowsize+1;
        
        endpoint=start+windowsize-1;
        
        if endpoint>N && length(start:endpoint)<ceil(windowsize/2)
            continue;
        elseif endpoint>N && length(start:endpoint)>= ceil(windowsize/2)
            temp=signal(:,start:end);
        else
            temp=signal(:,start:endpoint);  
        end
        
        
        temp=bsxfun(@minus, temp, mean(temp,2));
        temp=bsxfun(@rdivide, temp, std(temp, [],2));
        z_scored_signal=[z_scored_signal temp];
    end
    
end
 
    