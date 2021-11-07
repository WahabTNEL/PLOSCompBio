function [estParams, seq, LL]=GPFAtrain(data,xDim)
    

  binWidth      = 250; % in msec
  startTau      = 100; % in msec
  startEps      = 1e-3;
  
    % ==================================
  % Initialize state model parameters
  % ==================================
  startParams.covType = 'rbf';
  % GP timescale
  % Assume binWidth is the time step size.
  startParams.gamma = (binWidth / startTau)^2 * ones(1, xDim);
  % GP noise variance
  startParams.eps   = startEps * ones(1, xDim);
      yAll             = [data.y];
      [faParams, faLL] = fastfa(yAll, xDim);
  
  startParams.d = mean(yAll, 2);
  startParams.C = faParams.L;
  startParams.R = diag(faParams.Ph);

  % Define parameter constraints
  startParams.notes.learnKernelParams = true;
  startParams.notes.learnGPNoise      = false;
  startParams.notes.RforceDiagonal    = true;

  currentParams = startParams;
  
  [estParams, seq, LL, iterTime] = em(currentParams, data);
  
  
end
