function [prediction,accuracy]=ClassifyFactors(seq_train, y_train, seq_test,y_test, tlength, f_idx, type)
    
    n_dim=length(f_idx);
    for i=1:length(seq_train)
        train_features(i,:)=reshape(seq_train(i).xsm(f_idx,:), [n_dim*tlength 1]);
    end

    for i=1:length(seq_test)
        test_features(i,:)=reshape(seq_test(i).xsm(f_idx,:), [n_dim*tlength 1]);
    end

   %% Contruct parametrized mu and sigma
    K=length(unique(y_train));
    
  
    if strcmp(type, 'quadratic')==1
        
        
        for k=1:K
            
            train_data=[seq_train(find(y_train==k)).xsm];
            mu_small=mean(train_data(f_idx,:),2);
            mu_big(:,k)=repmat(mu_small, [tlength 1]);
            cov_small=cov(train_data(f_idx,:)');
            
            for t=1:tlength
                
                cov_big(1+(t-1)*n_dim:t*n_dim,1+(t-1)*n_dim:t*n_dim,k)=cov_small;
                
            end
        end
    
        Model = makecdiscr(mu_big',cov_big);
        
        
    elseif strcmp(type, 'linear')==1
        
         for k=1:K
             
            train_data=[seq_train(find(y_train==k)).xsm];
            mu_small=mean(train_data(f_idx,:),2);
            mu_big(:,k)=repmat(mu_small, [tlength 1]);
            
         end
         
            train_data=[seq_train.xsm];
            cov_small=cov(train_data(f_idx,:)');
            
        for t=1:tlength
            
            cov_big(1+(t-1)*n_dim:t*n_dim,1+(t-1)*n_dim:t*n_dim)=cov_small;
            
        end
    
        Model = makecdiscr(mu_big', cov_big);
    end
    
    [prediction,score] = predict(Model,test_features);
    numofcorrect=length(find((prediction-y_test')==0));


accuracy=numofcorrect./length(prediction);
end

