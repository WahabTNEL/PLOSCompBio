function method = plotPredErrorVsDim_v2(runIdxs, kernSD, varargin)
%
% method = plotPredErrorVsDim(runIdx, kernSD,...)
%
% Plot prediction error versus state dimensionality.
%
% INPUTS:
%
% runIdx    - results files will be loaded from mat_results/runXXX, where
%             XXX is runIdx
% kernSD    - smoothing kernel standard deviation to use for two-stage methods
%
% OUTPUTS:
%
% method    - data structure containing prediction error values shown in plot
%
% OPTIONAL ARGUMENTS:
%
% plotOn    - logical that specifies whether or not to display plot
%
% @ 2009 Byron Yu -- byronyu@stanford.edu

for R=1:length(runIdxs)
    runIdx=runIdxs(R);
  plotOn = true;
  assignopts(who, varargin);

  allMethods = {'pca', 'ppca', 'fa', 'gpfa'};
  
  runDir = sprintf('mat_results/run%03d', runIdx);
  if ~isdir(runDir)
    fprintf('ERROR: %s does not exist.  Exiting...\n', runDir);
    return
  else    
    D = dir([runDir '/*.mat']);
  end

  if isempty(D)
    fprintf('ERROR: No valid files.  Exiting...\n');
    method = [];
    return;
  end

  for i = 1:length(D)
    P = parseFilename(D(i).name);
    
    D(i).method = P.method;
    D(i).xDim   = P.xDim;
    D(i).cvf    = P.cvf;
    
    [tf, D(i).methodIdx] = ismember(D(i).method, allMethods);
  end
  % Only continue processing files that have test trials
  D = D([D.cvf]>0);

  if isempty(D)
    fprintf('ERROR: No valid files.  Exiting...\n');
    method = [];
    return;
  end
  
  for i = 1:length(D)    
    fprintf('Loading %s/%s...\n', runDir, D(i).name);
    ws = load(sprintf('%s/%s', runDir, D(i).name));
    
    % Check if selected kernSD has been run.
    if isfield(ws, 'kernSDList')
      % PCA, PPCA, FA
      [kernSD_exists, kidx] = ismember(kernSD, ws.kernSDList);
    else
      % GPFA
      kernSD_exists = true;
    end
    
    D(i).isValid = false;
    
    if kernSD_exists
      D(i).isValid = true;
      
      % Compute prediction error
      if D(i).methodIdx <= 3
        % PCA, PPCA, FA
        
        for c=1:max(ws.Y_test)
            
            Ycs            = [ws.kern(kidx).seqTest(find(ws.Y_test==c)).ycs];
            YtestRaw       = [ws.kern(kidx).seqTest(find(ws.Y_test==c)).y];
            %D(i).sse       = sum((Ycs(:)-ws.YtestRaw(:)).^2);
            T=length(Ycs(:));
            sse(c)=sqrt(sum((Ycs(:)-YtestRaw(:)).^2)/T);
            
        end
        
        D(i).sse=sum(sse);
        D(i).numTrials = length(ws.kern(kidx).seqTest);
      elseif D(i).methodIdx == 4
          
       
        % GPFA
        
        for p = 1:D(i).xDim
            
            for c=1:max(ws.Y_test)
                fn              = sprintf('ycsOrth%02d', p);
                Ycs             = [ws.seqTest(find(ws.Y_test==c)).(fn)];
                YtestRaw       = [ws.seqTest(find(ws.Y_test==c)).y];
            %  D(i).sseOrth(p) = sum((Ycs(:) - YtestRaw(:)).^2); 
                T=length(Ycs(:));
                sseOrth_c(c)=sqrt(sum((Ycs(:)-YtestRaw(:)).^2)/T);
                
            end
            
            D(i).sseOrth(p)=sum(sseOrth_c);
        end
        D(i).sse       = D(i).sseOrth(end);
        D(i).numTrials = length(ws.seqTest);        
      end
    end
  end
  
  D = D([D.isValid]);

  if isempty(D)
    fprintf('ERROR: No valid files.  Exiting...\n');
    method = [];
    return;
  end
  
  % Sum prediction error across cross-validation folds
  for n = 1:4
    Dn = D([D.methodIdx]==n);

    method(n).name = upper(allMethods{n});
    method(n).xDim = unique([Dn.xDim]);
    
    % Do for each unique state dimensionality.
    % Each method is allowed to have a different list of 
    % unique state dimensionalities.
    for p = 1:length(method(n).xDim)
      Dnn = Dn([Dn.xDim] == method(n).xDim(p));
      
      % Sum across cross-validation folds
      method(n).sse(p)       = mean([Dnn.sse]); 
      method(n).sse_sem(p)      = std([Dnn.sse])/sqrt(length([Dnn.sse]));
      method(n).numTrials(p) = sum([Dnn.numTrials]);
    end
  end
  
  % Reduced GPFA (based on GPFA files with largest xDim)
  Dn = D([D.methodIdx]==4);

  if length(Dn)>0
    dList = [Dn.xDim];
    Dnn   = Dn(dList == max(dList));
    
    method(5).name      = 'Reduced GPFA';
    method(5).xDim      = 1:max(dList);
    method(5).sse       = mean(vertcat(Dnn.sseOrth));
    method(5).sse_sem = std(vertcat(Dnn.sseOrth))/sqrt(size(vertcat(Dnn.sseOrth),1));
    method(5).numTrials = sum([Dnn.numTrials]);
  end
    
   FA(R,:)=method(3).sse;
  GPFA(R,:)=method(4).sse;
  RGPFA(R,:)=method(5).sse;
  numTrialsAll = [method.numTrials];
end
xDim=method(3).xDim;
 % if length(unique(numTrialsAll)) ~= 1
 %   fprintf('ERROR: Number of test trials must be the same across\n');
 %   fprintf('all methods and state dimensionalities.  Exiting...\n');
 %   return
 % end
    
  % =========
  % Plotting
  % =========
 %%
    col = {'r--', 'r', 'g', 'k--', 'k'}; 
    
    figure;

        plot(xDim, (FA), col{3}, 'LineWidth', 1.5);
        hold on;
                plot(xDim, (GPFA), col{4}, 'LineWidth', 1.5);
        hold on;
                plot(1:max(xDim), (RGPFA), col{5}, 'LineWidth', 1.5);
        hold on;
        %shadedErrorBar(method(n).xDim, method(n).sse, method(n).sse_sem, 'lineProps', col{n})

    legend('FA', 'GPFA', 'Reduced GPFA');
    legend('boxoff');
    %title(sprintf('For two-stage methods, kernel width = %d ms', kernSD),...
    %'fontsize', 12);
    xlabel('State dimensionality', 'fontsize', 14);
    ylabel('Prediction Error', 'fontsize', 14);
    xlim([1 max(xDim)]);
  end