function [Data_new] =RemoveZeros(Data)
    
N=length(Data);
T=size(Data{1},2);
C=size(Data{1},1);
    for i=1:N
         
        data_cell=Data{i};
        data_sum=sum(data_cell);
        
        zero_locs=find(data_sum==0);
        
        %Remove zeros
        
        data_cell(:,zero_locs)=[];
        
        Data_new{i}=data_cell;
    end
end

            
