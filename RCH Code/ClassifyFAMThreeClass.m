function [C, accuracy]=ClassifyFAMThreeClass(DTrainhil, RTrainhil, ELTrainhil,  DTesthil, RTesthil, ELTesthil, numofFactors, ifFactorAnalysis)



DPos=DTrainhil;
RPos=RTrainhil;
ELPos=ELTrainhil;



TrainingSet=[DPos RPos ELPos]';
LabelTr=[1*ones(1,size(DPos,2)) 2*ones(1,size(RPos,2)) 3*ones(1,size(ELPos,2))];

m1=mean(DPos,2);
m2=mean(RPos,2);
m3=mean(ELPos,2);

n1=length(DPos);
n2=length(RPos);
n3=length(ELPos);


s1=var(DPos, [], 2);
s2=var(RPos, [], 2);
s3=var(ELPos, [], 2);



mt=mean(TrainingSet,1)';
nt=length(TrainingSet);

meanTrain=(m1*n1+m2*n2+m3*n3)/nt;
stdTrain=sqrt((n1*s1+n2*s2+n3*s3+n1*(m1-mt).^2+n2*(m2-mt).^2+n3*(m3-mt).^2)/nt);

%[~, idx]=licols(TrainingSet);
meanTrain=mean(TrainingSet, 1);
stdTrain=sqrt(var(TrainingSet, [], 1));
%TrainingSet=bsxfun(@minus, TrainingSet, meanTrain);
%TrainingSet=bsxfun(@rdivide, TrainingSet, stdTrain);
%TrainingSet=TrainingSet(:,idx);
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



TestSet=[DPosTs RPosTs ELPosTs]';

LabelTs=[ones(1,size(DPosTs,2)) 2*ones(1,size(RPosTs,2)) 3*ones(1,size(ELPosTs,2))];

%TestSet=bsxfun(@minus, TestSet, mean(TestSet, 1));
%TestSet=bsxfun(@rdivide, TestSet, sqrt(var(TestSet, [], 1)));

%TestSet=bsxfun(@minus, TestSet, meanTrain);
%TestSet=bsxfun(@rdivide, TestSet, stdTrain);

%TestSet=TestSet(:,idx);
t=templateSVM('KernelFunction', 'linear');
Model=fitcecoc(mappedX', LabelTr,  'Learners', t);



if strcmp(ifFactorAnalysis, 'Yes')
mappedXTs=lambda'*inv(lambda*lambda'+Psi)*TestSet';
else
mappedXTs=TestSet';
end
[prediction,score] = predict(Model,mappedXTs');

numofcorrect=length(find((prediction-LabelTs')==0));


accuracy=numofcorrect./length(prediction);

[C, order] = confusionmat(LabelTs, prediction);
C(1,:) = C(1,:)./sum(C(1,:));
C(2,:) = C(2,:)./sum(C(2,:));
C(3,:) = C(3,:)./sum(C(3,:));



end