function result = neuralTraj_RCH_v2(runIdx, dat, indD, indR, indEL,  varargin)
%
% result = neuralTraj(runIdx, dat, ...)
%
% Prepares data and calls functions for extracting neural trajectories.
%
% INPUTS:
%
% runIdx      - results files will be saved in mat_results/runXXX, where
%               XXX is runIdx
% dat         - structure whose nth entry (corresponding to the nth experimental
%               trial) has fields
%                 trialId -- unique trial identifier
%                 spikes  -- 0/1 matrix of the raw spiking activity across 
%                            all neurons.  Each row corresponds to a neuron.  
%                            Each column corresponds to a 1 msec timestep.
%
% OUTPUTS:
%
% result      - structure containing all variables saved in mat_results/runXXX/
%               if 'numFolds' is 0.  Else, the structure is empty.
%               
% OPTIONAL ARGUMENTS:
%
% method      - method for extracting neural trajectories
%               'gpfa' (default), 'fa', 'ppca', 'pca'
% binWidth    - spike bin width in msec (default: 20)
% numFolds    - number of cross-validation folds (default: 0)
%               0 indicates no cross-validation, i.e. train on all trials.
% xDim        - state dimensionality (default: 3)
%
% @ 2009 Byron Yu         byronyu@stanford.edu
%        John Cunningham  jcunnin@stanford.edu

  method        = 'gpfa'; 
  binWidth      = 250; % in msec
  numFolds      = 0;
  xDim          = 8;
  K=[];
  Y=[];
  extraOpts     = assignopts(who, varargin);

  fprintf('\n---------------------------------------\n');
  if ~isdir('mat_results')
    mkdir('mat_results');
  end
  % Make a directory for this runIdx if it doesn't already exist
  runDir = sprintf('mat_results/run%03d', runIdx);
  if isdir(runDir)
    fprintf('Using existing directory %s...\n', runDir);
  else
    fprintf('Making directory %s...\n', runDir);
    mkdir(runDir);
  end

  % Obtain binned spike counts
%  seq  = getSeq(dat, binWidth, extraOpts{:});
seq=dat;
  if isempty(seq)
    fprintf('Error: No valid trials.  Exiting.\n');
    result = [];
    return;
  end
  
  % Set cross-validation folds 
  numFolds=K;

  
  seq_d=seq(find(Y==0));
  seq_r=seq(find(Y==1));
  seq_el=seq(find(Y==2));
  
  
  for cvf = 0:numFolds
      
    if cvf == 0      
      fprintf('\n===== Training on all data =====\n');
    else
      fprintf('\n===== Cross-validation fold %d of %d =====\n', cvf, numFolds);
    end

    % Specify filename where results will be saved
    fname = sprintf('%s/%s_xDim%02d', runDir, method, xDim);
    if cvf > 0
      fname = sprintf('%s_cv%02d', fname, cvf);
    end
    if exist([fname '.mat'], 'file')
      fprintf('%s already exists.  Skipping...\n', fname);
      continue;
    end
    
    % Set cross-validation train and test sets
    d=1;

    if cvf > 0
        
        
          seq_d_train1=seq_d(find(indD<cvf-d | indD>cvf+d));
          seq_r_train1=seq_r(find(indR<cvf-d | indR>cvf+d));
          seq_el_train1=seq_el(find(indEL<cvf-d | indEL>cvf+d));
          
          minlength=min([length(seq_d_train1) length(seq_r_train1) length(seq_el_train1)]);
          seq_d_train=seq_d_train1(randperm(length(seq_d_train1), minlength));
          seq_r_train=seq_r_train1(randperm(length(seq_r_train1), minlength));
          seq_el_train=seq_el_train1(randperm(length(seq_el_train1), minlength));
          
          seqTrain=[seq_d_train seq_r_train seq_el_train];
          
          seq_d_test1=seq_d(find(indD==cvf));
          seq_r_test1=seq_r(find(indR==cvf));
          seq_el_test1=seq_el(find(indEL==cvf));
          
          
          seqTest=[seq_d_test1 seq_r_test1 seq_el_test1];          
          Y_test=[ones(1,length(seq_d_test1)) 2*ones(1, length(seq_r_test1)) 3*ones(1,length(seq_el_test1))];
          
    else
          seq_d_train1=seq_d;
          seq_r_train1=seq_r;
          seq_el_train1=seq_el;
          
          minlength=min([length(seq_d_train1) length(seq_r_train1) length(seq_el_train1)]);
          seq_d_train=seq_d_train1(randperm(length(seq_d_train1), minlength));
          seq_r_train=seq_r_train1(randperm(length(seq_r_train1), minlength));
          seq_el_train=seq_el_train1(randperm(length(seq_el_train1), minlength));
          
          seqTrain=[seq_d_train seq_r_train seq_el_train];
          %Create dummy empty struct for test data
          seqTest=seq(find(Y==Inf));
          Y_test=[];
    end
    

      
          % If training on all trials, keep original trial ordering

        
        %trainTrialIdx = tr(trainMask);
        %testTrialIdx  = tr(testMask);    
        %seqTrain      = seq(trainTrialIdx);
        %seqTest       = seq(testTrialIdx);
    
    
 
    
    % Remove inactive units based on training set
    hasSpikesBool = (mean([seqTrain.y], 2) ~= 0);
  
    for n = 1:length(seqTrain)
      seqTrain(n).y = seqTrain(n).y(hasSpikesBool,:);
    end
    
    
    for n = 1:length(seqTest)
      seqTest(n).y = seqTest(n).y(hasSpikesBool,:);
    end

    % Check if training data covariance is full rank
    yAll = [seqTrain.y];
    yDim  = size(yAll, 1);
    
    if rank(cov(yAll')) < yDim
      fprintf('ERROR: Observation covariance matrix is rank deficient.\n');
      fprintf('Possible causes: repeated units, not enough observations.\n');
      fprintf('Exiting...\n');
      return
    end
    
    fprintf('Number of training trials: %d\n', length(seqTrain));
    fprintf('Number of test trials: %d\n', length(seqTest));
    fprintf('Latent space dimensionality: %d\n', xDim);
    fprintf('Observation dimensionality: %d\n', sum(hasSpikesBool));

    % If doing cross-validation, don't use private noise variance floor.  
    if cvf > 0
      extraOpts = {extraOpts{:}, 'minVarFrac', -Inf};      
    end

    % The following does the heavy lifting.
    if isequal(method, 'gpfa')
      gpfaEngine_v2(seqTrain, seqTest, Y_test, fname,... 
      'xDim', xDim, 'binWidth', binWidth, extraOpts{:});
    
    elseif ismember(method, {'fa', 'ppca', 'pca'})
      twoStageEngine_v2(seqTrain, seqTest, Y_test, fname,... 
      'typ', method, 'xDim', xDim, 'binWidth', binWidth, extraOpts{:});
    end

    if exist([fname '.mat'], 'file')
      save(fname, 'method', 'cvf', 'hasSpikesBool', 'extraOpts', '-append');
    end
  end
  
  result = [];  
  if (nargout == 1) & (numFolds == 0) & exist([fname '.mat'], 'file')
    result = load(fname);
  end
