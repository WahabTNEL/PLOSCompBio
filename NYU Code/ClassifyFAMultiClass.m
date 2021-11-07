function [C, accuracy]=ClassifyFAMultiClass(DTrainhil, RTrainhil, ELTrainhil, TVTrainhil,  DTesthil, RTesthil, ELTesthil, TVTesthil, numofFactors, ifFactorAnalysis)

DPos=DTrainhil;
RPos=RTrainhil;
ELPos=ELTrainhil;
TVPos=TVTrainhil;


TrainingSet=[DPos RPos ELPos TVPos]';

LabelTr=[ones(1,size(DPos,2)) 2*ones(1,size(RPos,2)) 3*ones(1,size(ELPos,2)) 4*ones(1, size(TVPos,2))];

%%

 if strcmp(ifFactorAnalysis, 'Yes')
%[TrainingSet, idx]=licols(TrainingSet);
[~, mapping, Psi]=fa(TrainingSet, numofFactors);
%[lambda, Psi]=factoran(TrainingSet, numofFactors);

%Psi=diag(Psi);
lambda=mapping.M;
mappedX=lambda'*inv(lambda*lambda'+Psi)*TrainingSet';
else
    mappedX=TrainingSet';
end
%%
DPosTs=DTesthil;
RPosTs=RTesthil;
ELPosTs=ELTesthil;
TVPosTs=TVTesthil;

    TestSet=[DPosTs RPosTs ELPosTs TVPosTs]';
    LabelTs=[ones(1,size(DPosTs,2)) 2*ones(1,size(RPosTs,2)) 3*ones(1,size(ELPosTs,2)) 4*ones(1, size(TVPosTs,2))];
  
    

Model=fitcecoc(mappedX', LabelTr);
%obj=gmdistribution.fit(mappedX', 4);

if strcmp(ifFactorAnalysis, 'Yes')
mappedXTs=lambda'*inv(lambda*lambda'+Psi)*TestSet';
%mappedXTs=mapping.M'*TestSet';
else
    mappedXTs=TestSet';
end


[prediction,score] = predict(Model,mappedXTs');
%prediction=cluster(obj, mappedXTs');
numofcorrect=length(find((prediction-LabelTs')==0));


accuracy=numofcorrect./length(prediction);

[C, order] = confusionmat(LabelTs, prediction);
C(1,:) = C(1,:)./sum(C(1,:));
C(2,:) = C(2,:)./sum(C(2,:));
C(3,:) = C(3,:)./sum(C(3,:));
C(4,:) = C(4,:)./sum(C(4,:));

%figure;
%plotCM(C, 0,{'Dialogue', 'Rest'})

end