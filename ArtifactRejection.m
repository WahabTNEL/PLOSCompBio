function [bad_epochs] =ArtifactRejection(Data, param1, param2, param3)

N=length(Data);
T=size(Data{1},2);
C=size(Data{1},1);
bad_epochs=zeros(1,N);
for e=1:N
    
      epoch=Data{e};
      ChannelMean=mean(epoch,2);
      ChannelStd=std(epoch, [], 2);
      bad_time=zeros(1,T);
         for t=1:T
             for c=1:C
                 bad_channel(c)=epoch(c,t)>(ChannelMean(c)+ChannelStd(c)*param1) || epoch(c,t)<(ChannelMean(c)-ChannelStd(c)*param1);
             end  
          if (sum(bad_channel)/length(bad_channel))>param2 
             bad_time(t)=1;             
          end
         end
         
      if sum(bad_time) >= param3
          bad_epochs(e)=1;
      end
      
end