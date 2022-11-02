classdef dataset
% DATASET Class containing a dataset (usually a list of experiments from a
% given animal)
%
% USAGE:
%   obj = dataset(datasetName, varargin);
%
% INPUT arguments:
%   datasetName - Name of the dataset to create or load (if it is an
%   existing .mat file)
%
% INPUT optional arguments ('key' followed by its value):
%   rootFolder - folder where to create/look for the dataset (current
%   directory by default
%   dataFolder - where to store data (rootFolder/datasetName/data by default)
%   outputFolder - where to store the output (rootFolder/datasetName/output by default)
%   comments - dataset comments you might want to add
%
% OUTPUT arguments:
%   obj - dataset structure
%
% EXAMPLE:
%   dataset = dataset('animal1', 'baseFolder', pwd);
%
% Copyright (C) 2019-2021, Javier G. Orlandi <javierorlandi@javierorlandi.com>,
% Benuccilab (http://benuccilab.brain.riken.jp)
%
% See also class methods: save
  
  properties
    name;
    rootFolder;
    dataFolder;
    outputFolder;
    comments;
    tags; % List of tags to filter/sort stuff
    info; % Struct to be used for any other subclass
    data; % Structure that contains data
  end
  methods(Static)
    function obj = load(rootFolder, datasetName, varargin)
      params.verbose = true;
      params = parse_pv_pairs(params, varargin);
      obj = wf.dataset();
      % Assign loaded properties to the class
      if(exist(fullfile(rootFolder, datasetName, [datasetName '.mat']), 'file'))
        if(params.verbose)
          slog('Loading: %s', datasetName);
        end
        data = load(fullfile(rootFolder, datasetName, [datasetName '.mat'])); %#ok<*PROPLC>
        fnames = fieldnames(data);
        for it = 1:length(fnames)
          try
            obj.(fnames{it}) = data.(fnames{it});
          catch ME
            slog.warning(ME.message);
          end
        end
      else
        slog.warning('Dataset: %s not found', datasetName);
      end
    end
  end
  methods
    % The actual constructor
    function obj = dataset(datasetName, varargin)
      if(nargin == 0)
        obj.name = 'empty';
        obj.rootFolder = [];
        obj.dataFolder = [];
        obj.outputFolder = [];
        obj.comments = 'empty dataset';
        obj.tags = {};
        obj.info = struct;
        obj.data = struct;
        return;
      end
      % DATASET class constructor
      params.rootFolder = pwd;
      params.dataFolder = 'data';
      params.outputFolder = 'output';
      params.comments = '';
      params.tags = '';
      params.continuous = false;
      params.verbose = true;
      params = parse_pv_pairs(params, varargin);
      
      % Check if we are passing an existing file
%       if(exist(datasetName, 'file') == 2)
%         [fpa, fpb, fpc] = fileparts(datasetName);
%         if(strcmpi(fpc, '.mat'))
%           slog.info('Dataset %s already exists', fpb);
%           obj = obj.load(params.baseFolder, datasetName);
%           return;
%         end
      % Check if we the target directory already exists and contains a
      % dataset
      if(exist(fullfile(params.rootFolder, datasetName, [datasetName '.mat']), 'file'))
        slog.info('Dataset %s already exists', datasetName);
        obj = obj.load(params.rootFolder, datasetName);
        return;
      elseif(exist(datasetName, 'file') == 2)
        [fpa, fpb, ~] = fileparts(datasetName);
        datasetName = fpb;
        params.rootFolder = string(unique(arrayfun(@(x) x.folder, dir([fpa filesep '..' filesep '..']), 'UniformOutput', false)));
        slog.info('Dataset %s already exists', datasetName);
        obj = obj.load(params.rootFolder, datasetName);
        return;
      end
      
      % Else, create everything
      obj.name = datasetName;
      obj.rootFolder = params.rootFolder;
      obj.dataFolder = params.dataFolder;
      obj.outputFolder = params.outputFolder;
      obj.comments = params.comments;
      obj.info.continuous = params.continuous;
      obj.tags = params.tags;
      if(params.verbose)
        slog.info('Creating dataset %s on: %s', datasetName, obj.rootFolder);
      end
      
      % Create the required folders
      baseFolder = fullfile(obj.rootFolder, datasetName);
      [~, msg] = mkdir(baseFolder);
      if(~isempty(msg))
        slog.warning('mkdir %s: %s', baseFolder, msg);
      end
      [~, msg] = mkdir(fullfile(baseFolder,obj.dataFolder));
      if(~isempty(msg))
        slog.warning('mkdir %s: %s', fullfile(baseFolder,obj.dataFolder), msg);
      end
      [~, msg] = mkdir(fullfile(baseFolder,obj.outputFolder));
      if(~isempty(msg))
        slog.warning('mkdir %s: %s', fullfile(baseFolder,obj.outputFolder), msg);
      end
      % Save the dataset
      obj.save();
    end

    function obj = loadTrialsSummary(obj, fileName, varargin)
    % LOADTRIALSSUMMARY Loads the trials summary information from file
    %
    % USAGE:
    %   obj = loadTrialsSummary(trialsSummaryFile);
    % INPUT arguments:
    %   trialsSummaryFile - file containing the trials summary information 
    %                       (the L file)
    %
    % INPUT optional arguments ('key' followed by its value):
    %   verbose - displays additional info (true/false)
    %
    % OUTPUT arguments:
    %   obj - class object
    %
    % EXAMPLE:
    %   dataset = wf.dataset.load(pwd, 'animal1');
    %   dataset = dataset.loadTrialsSummary('S.mat');
    
      params.verbose = true;
      params = parse_pv_pairs(params, varargin);
      
      tmp = load(fileName);
      obj.data.trialsSummary = tmp.L;
      [fpa, ~, ~] = fileparts(fileName);
      obj.info.trialsSummaryFolder = fpa;
      if(isfield(obj.info, 'validSessions') && isempty(obj.info.validSessions))
        obj.info.validSessions = unique(obj.data.trialsSummary.sess);
      end
      if(params.verbose)
        slog.info('Trials summary data succesfully loaded');
      end
    end
    
    function obj = alignToReference(obj, referenceDataset, varargin)
      params.verbose = true;
      params = parse_pv_pairs(params, varargin);
      
      ref = wf.dataset.load(obj.rootFolder, referenceDataset, 'verbose', false);
      refTriangle = ref.data.retinotopyCoordinates.alignmentTriangle; % The first of them
      refName = ref.name;

      alignmentTriangle = obj.data.retinotopyCoordinates.alignmentTriangle;
  
      trTV1ref = eye(3); % Translate to global ref orientation
      trTC = eye(3); % Translate to centroid
      trS = eye(3); % Scale from centroid
      trTMC = eye(3); % Translate back
      trTV1 = eye(3); % Translate to V1
      trR = eye(3); % Rotate around V1
      trTMV1 = eye(3); % Translate back

      trFull = affine2d;

      % Translate V1 to reference V1
      trTV1ref(3, [1 2]) = refTriangle.Vertices(1,:) - alignmentTriangle.Vertices(1, :);
      % Translate V1 to 0
      trTC(3, [1 2]) = -refTriangle.Vertices(1,:);
      % Scale to match reference perimeter
      trS(1, 1) = perimeter(refTriangle)/perimeter(alignmentTriangle);
      trS(2, 2) = perimeter(refTriangle)/perimeter(alignmentTriangle);
      % Go back to reference V1
      trTMC(3, [1 2]) = refTriangle.Vertices(1,:);

      % move the new V1 to 0 again for the rotation
      trFull.T = trTV1ref*trTC*trS*trTMC;
      [xn,yn] = transformPointsForward(trFull,alignmentTriangle.Vertices(:,1), alignmentTriangle.Vertices(:,2));
      trTV1(3, [1 2]) = -[xn(1) yn(1)];
      trTMV1(3, [1 2]) = [xn(1) yn(1)];
      % move V1 of the reference triangle to 0
      nRefTriangle = refTriangle.Vertices-[xn(1) yn(1)];
      % Need to know where the triangle points are after centering at V1
      trFull.T = trTV1ref*trTC*trS*trTMC*trTV1;
      [xn,yn] = transformPointsForward(trFull,alignmentTriangle.Vertices(:,1), alignmentTriangle.Vertices(:,2));

      % Move to polar coordinates
      xr = nRefTriangle(2:3, 1);
      yr = nRefTriangle(2:3, 2);
      m1 = norm([xn(2), yn(2)]);
      m2 = norm([xn(3), yn(3)]);
      t1 = atan2(yn(2),xn(2));
      t2 = atan2(yn(3),xn(3));
      % Distance minimizer function
      mfun = @(x) norm([xr(1)-m1*cos(t1+x) yr(1)-m1*sin(t1+x)])+ norm([xr(2)-m2*cos(t2+x), yr(2)-m2*sin(t2+x)]);
      [mt, fv] = fminbnd(mfun, -pi/4, pi/4);

      % Compute the rotation
      trR(1:2, 1:2) = [cos(mt) sin(mt);-sin(mt) cos(mt)];
      % Do the full transformation - after the rotation I have to move V1 back to its position
      trFull.T = trTV1ref*trTC*trS*trTMC*trTV1*trR*trTMV1; % Need to know where the triangle points are after centering at V1
      [xn,yn] = transformPointsForward(trFull,alignmentTriangle.Vertices(:,1), alignmentTriangle.Vertices(:,2));

      obj.data.alignment = struct;
      obj.data.alignment.transform = trFull;
      obj.data.alignment.err = fv;
      obj.data.alignment.referenceName = refName;
      obj.data.alignment.referenceTriangle = refTriangle;
      obj.data.alignment.triangle = polyshape(xn', yn');
      
      if(params.verbose)
        slog.info('Alignment to reference completed.');
      end
    end
    
    function sessList = getValidSessions(obj)
    % GETVALIDSESSIONS Gets the list of sessions that will be used during
    % preprocessing
    %
    % USAGE:
    %   sessList = obj.getValidSessions();
    %
    % INPUT arguments:
    %   none
    %
    % OUTPUT arguments:
    %   sessList - vector containing the session IDs that will be used
    %
    % EXAMPLE:
    %   dset = wf.dataset.load(pwd, 'animal1');
    %   sessList = dset.getValidSessions();
      
      sessList = obj.info.validSessions(:)';
    end
    
    function obj = setValidSessions(obj, sessList)
    % SETVALIDSESSIONS Sets the list of sessions that will be used during
    % preprocessing
    %
    % USAGE:
    %   obj = obj.setValidSessions(sess);
    %
    % INPUT arguments:
    %   sessList - vector containing the session IDs to use (subset of the
    %              obj.data.trialSummary.sess list)
    %
    % OUTPUT arguments:
    %   obj - class object
    %
    % EXAMPLE:
    %   dset = wf.dataset.load(pwd, 'animal1');
    %   dset = dset.setValidSessions(1:5);
      
      obj.info.validSessions = sessList(:);
    end
    
    function obj = setSessionList(obj, sessList)
    % SETSESSIONLIST Sets the session list IDs
    %
    % USAGE:
    %   obj = obj.setSessionList(sess);
    %
    % INPUT arguments:
    %   sessList - cell containing the sessions IDS {'name', 'sess', 'exp'}
    %
    % OUTPUT arguments:
    %   obj - class object
    %
    % EXAMPLE:
    %   dset = wf.dataset.load(pwd, 'animal1');
    %   dset = dset.setValidSessions(1:5);
      
      obj.info.sessionList = sessList;
    end
    
    function obj = resetValidSessions(obj)
    % RESETVALIDSESSIONS Resets the list of sessions that will be used
    % during preprocessing (from the trialsSummary file)
    %
    % USAGE:
    %   obj = obj.resetValidSessions();
    %
    % INPUT arguments:
    %   none
    %
    % OUTPUT arguments:
    %   obj - class object
    %
    % EXAMPLE:
    %   dset = wf.dataset.load(pwd, 'animal1');
    %   dset = dset.resetValidSessions();
      if(isfield(obj.data, 'trialsSummary'))
        obj.info.validSessions = unique(obj.data.trialsSummary.sess);
      else
        obj.info.validSessions = [];
      end
    end
    
    function success = save(obj)
    % SAVE Save the dataset to file
    %
    % USAGE:
    %   success = obj.save();
    %
    % INPUT arguments:
    %   none
    %
    % OUTPUT arguments:
    %   success - whether the save was succesfull (true/false)
    %
    % EXAMPLE:
    %   dset = wf.dataset.load(pwd, 'animal1');
    %   success = dset.save();
    
      s = struct;
      fnames = fieldnames(obj);
      for it = 1:length(fnames)
         try
          s.(fnames{it}) = obj.(fnames{it});
        catch ME
          slog.warning(ME.message);
         end
      end
      slog('Saving: %s', fullfile(obj.rootFolder, obj.name, [obj.name '.mat']));
      try
        save(fullfile(obj.rootFolder, obj.name, [obj.name '.mat']), '-struct', 's');
        success = true;
      catch
        success = false;
      end
    end
    
    %% --------------------------------------------------------------------
    %%% Plotting methods
    %%% -------------------------------------------------------------------
    
    function exportFig(obj, hFig, params)
      ext = params.export;
      exportFile = sprintf('%s%s_%s.%s', params.mainTag, params.appendTag, obj.name, ext);
      exportFileSub = sprintf('%s_%s%s.%s', obj.name, params.mainTag, params.appendTag, ext);
      if(~exist(fullfile(obj.rootFolder, obj.outputFolder, 'plots'), 'dir'))
        mkdir(fullfile(obj.rootFolder, obj.outputFolder, 'plots'))
      end
      if(~exist(fullfile(obj.rootFolder, obj.outputFolder, 'plots', obj.name), 'dir'))
        mkdir(fullfile(obj.rootFolder, obj.outputFolder, 'plots', obj.name))
      end
      exportgraphics(hFig, fullfile(obj.rootFolder, obj.outputFolder, 'plots', exportFile));
      copyfile(fullfile(obj.rootFolder, obj.outputFolder, 'plots', exportFile), ...
               fullfile(obj.rootFolder, obj.outputFolder, 'plots', obj.name, exportFileSub), 'f');
      if(params.verbose)
        slog('%s succesfully generated in %s', exportFile, fullfile(obj.rootFolder, obj.outputFolder, 'plots'));
      end
    end
    %%
    function hFig = plotFps(obj, varargin)
    % PLOTFPS Plots the FPS of each trial (from the validSessions set)
    % USAGE:
    %   obj.plotFps();
    %
    % INPUT arguments:
    %   none
    %
    % INPUT optional arguments ('key' followed by its value):
    %   rounding  - decimal number to round to, e.g., 0 for no decimals,
    %               1, for the first decimal, etc. Empty for no rounding
    %               (default). ([], 0, 1, ...)
    %   export    - to automatically export the image to a pdf or png file.
    %               Empty for no export (default). ([], 'pdf', 'png')
    %   mainTag   - main tag to use for the export file name. Default: '_fps' (text)
    %   appendTag - tag to append to the end of the exported file name. Empty
    %               by default. ([], text)
    %   visible   - if the figure should be visible. Only useful for exports.
    %               If not visible, it will be closed  at the end of the 
    %               script (on by default). (true/false)
    %   verbose   - displays additional info (true/false)
    %
    % OUTPUT arguments:
    %   hFig - Handle of the new figure
    %
    % EXAMPLE:
    %   dset = wf.dataset.load(pwd, 'animal1');
    %   dset.plotFps('rounding', 3, 'export', 'pdf');
      
      params.rounding = [];
      params.export = [];
      params.appendTag = '';
      params.mainTag = 'fps';
      params.visible = true;
      params.verbose = true;
      params = parse_pv_pairs(params, varargin);
      if(params.visible)
        visible = 'on';
      else
        visible = 'off';
      end
      
      if(params.verbose)
        slog('Plotting fps: %s', obj.name);
      end
      
      sessList = obj.getValidSessions();

      validTrials = find(sum(obj.data.trialsSummary.sess == sessList,2));
      
      fps = obj.data.trialsSummary.fps(validTrials,1);
      if(~isempty(params.rounding))
        fps = round(fps*10^params.rounding)/10^params.rounding;
      end
      hFig = figure('Visible', visible, 'Color', 'w');
      plot(validTrials, fps);
      xlabel('Trial index');
      ylabel('fps');
      xrange = minmax(validTrials(:)');
      yrange = minmax(fps(:)');
      if(diff(xrange))
        xl = xrange+[-1 1]*diff(xrange)*0.025;
      else
        xl = xlim;
      end
      if(diff(yrange))
        yl = yrange+[-1 1]*diff(yrange)*0.025;
      else
        yl = ylim;
      end
      ax1 = gca;
      ax1_pos = ax1.Position; % position of first axes
      sessFirstTrial = arrayfun(@(x)find(obj.data.trialsSummary.sess(validTrials) == x, 1, 'first'), sessList');
      hold on;
      box off;
      plot([sessFirstTrial sessFirstTrial], yl, 'k:');
      ax2 = axes('Position',ax1_pos,...
                  'XAxisLocation','top',...
                  'YAxisLocation','right',...
                  'Color','none');
      ax2.YTick = [];
      ax2.XTick = sessFirstTrial;
      ax2.XTickLabels = obj.info.validSessions;
      xlabel('Session start');
      ax2.XLim = xl;
      ax2.YLim = yl;
      linkaxes([ax1 ax2]);
      title(sprintf('FPS across sessions/trials for: %s', obj.name));
      TightInset = max([ax1.TightInset; ax2.TightInset]); %[Left, Bot, Right, Top] border needed
      AxPos = [TightInset(1), TightInset(2), 1-sum(TightInset([1 3])), 1-sum(TightInset([2 4]))]; %New position
      ax1.Position = AxPos;
      ax2.Position = AxPos;
      ax2.XLim = xl;
      ax2.YLim = yl;
      ax1.XLim = xl;
      ax2.YLim = yl;
      
      if(~isempty(params.export))
        obj.exportFig(hFig, params);
      end
      if(~params.visible)
        close(hFig);
      end
    end
    %%
    function [hFig, data] = plotSessionPerformance(obj, varargin)
    % PLOTSESSIONPERFORMANCE Plots the performance for each session
    % USAGE:
    %   obj.plotSessionPerformance();
    %
    % INPUT arguments:
    %   none
    %
    % INPUT optional arguments ('key' followed by its value):
    %   doBarPlot - Plot the distribution of C/I/TO as a bar plot instead.
    %               (true/false (default))
    %   export    - to automatically export the image to a pdf or png file.
    %               Empty for no export. ([] (default), 'pdf', 'png')
    %   mainTag   - main tag to use for the export file name. (txt, '_fps' (default))
    %   appendTag - tag to append to the end of the exported file name. Empty
    %               by default. ([] (default), text)
    %   visible   - if the figure should be visible. Only useful for exports.
    %               If not visible, it will be closed  at the end of the 
    %               script (on by default). ('on','off')
    %   verbose   - displays additional info (true/false)
    %
    % OUTPUT arguments:
    %   hFig - Handle of the new figure
    %   data - Data from the curves (as a struct)
    %
    % EXAMPLE:
    %   dset = wf.dataset.load(pwd, 'animal1');
    %   dset.plotSessionPerformance('export', 'pdf');

      params.doBarPlot = false;
      params.export = [];
      params.appendTag = '';
      params.mainTag = 'sessperf';
      params.visible = true;
      params.verbose = true;
      params = parse_pv_pairs(params, varargin);
      if(params.visible)
        visible = 'on';
      else
        visible = 'off';
      end
      
      if(params.verbose)
        slog('Plotting session performance: %s', obj.name);
      end
      
      sessList = obj.getValidSessions();

      sessPerf = zeros(length(sessList), 1);
      sessPerfNoTO = zeros(length(sessList), 1);
      sessCITO = zeros(length(sessList), 3);
      sessTO = zeros(length(sessList), 1);
      for it0 = 1:length(sessList)
        
        sessPerf(it0) = sum(obj.data.trialsSummary.outcome(ismember(obj.data.trialsSummary.sess, sessList(it0))) == 1) / sum(ismember(obj.data.trialsSummary.sess, sessList(it0)));
        sessPerfNoTO(it0) = sum(obj.data.trialsSummary.outcome(ismember(obj.data.trialsSummary.sess, sessList(it0))) == 1) / sum(obj.data.trialsSummary.outcome(ismember(obj.data.trialsSummary.sess, sessList(it0))) == 1 | obj.data.trialsSummary.outcome(ismember(obj.data.trialsSummary.sess, sessList(it0))) == -1);
        sessTO(it0) = sum(obj.data.trialsSummary.outcome(ismember(obj.data.trialsSummary.sess, sessList(it0))) == 0) / sum(ismember(obj.data.trialsSummary.sess, sessList(it0)));
        sessCITO(it0, :) = [sum(obj.data.trialsSummary.outcome(ismember(obj.data.trialsSummary.sess, sessList(it0))) ==  1), ...
                            sum(obj.data.trialsSummary.outcome(ismember(obj.data.trialsSummary.sess, sessList(it0))) == -1), ...
                            sum(obj.data.trialsSummary.outcome(ismember(obj.data.trialsSummary.sess, sessList(it0))) ==  0)]/ sum(ismember(obj.data.trialsSummary.sess, sessList(it0)));
      end

      hFig = figure('Visible', visible, 'Color', 'w');
      hold on;
      if(params.doBarPlot)
        bar(sessList, sessCITO, 0.8, 'stacked');
        ylim([0 1]);
        xlim([min(sessList) max(sessList)]+[-0.4 0.4]);
        l = legend('C', 'I', 'TO');
        ylabel('outcome fraction');
      else
        plot(sessList, sessPerf);
        plot(sessList, sessPerfNoTO);
        plot(sessList, sessTO);
        xlim([min(sessList) max(sessList)]);
        ylim([0 1]);
        spaceOutAxes();
        l = legend('Performance', 'Perf excl TO', 'Time out rate');
        ylabel('Performance');
      end
      xlabel('Session');
      box on;

      if(params.doBarPlot)
        l.Orientation = 'horizontal';
        l.Location = 'northOutside';
      else
        l.Location = 'best';
      end
      l.Box = 'off';
      title(sprintf('Session performance and time out rate: %s', obj.name));
      if(~params.doBarPlot)
        spaceOutAxes();
        offsetAxes();
      end
      data = struct;
      data.sessList = sessList;
      data.perf = sessPerf;
      data.perfNoTO = sessPerfNoTO;
      data.sessTO = sessTO;
      data.sessCITO = sessCITO;
      if(~isempty(params.export))
        obj.exportFig(hFig, params);
      end
      if(~params.visible)
        close(hFig);
      end
    end
    %%
    function hFig = plotPsychometricCurves(obj, varargin)
    % PLOTPSYCHOMETRICCURVES Plots the psychometric curves for each session and the grand average
    % USAGE:
    %   obj.plotSessionPsychometricCurve();
    %
    % INPUT arguments:
    %   none
    %
    % INPUT optional arguments ('key' followed by its value):
    %   showIndividualSessions - If it should plot each individual session (default true)
    %   lineAlpha - Alpha level for the session lines (default 0.5)
    %   sessionList - Subset of sessions to use instead (default all valid sessions)
    %   angleDivisions - Number of divisions for the angle list
    %   export    - to automatically export the image to a pdf or png file.
    %               Empty for no export. ([] (default), 'pdf', 'png')
    %   mainTag   - main tag to use for the export file name. (txt, '_fps' (default))
    %   appendTag - tag to append to the end of the exported file name. Empty
    %               by default. ([] (default), text)
    %   visible   - if the figure should be visible. Only useful for exports.
    %               If not visible, it will be closed  at the end of the 
    %               script (on by default). (true/false)
    %   verbose   - displays additional info (true/false)
    %
    % OUTPUT arguments:
    %   hFig - Handle of the new figure
    %
    % EXAMPLE:
    %   dset = wf.dataset.load(pwd, 'animal1');
    %   dset.plotSessionPsychometricCurve('export', 'pdf');

      params.export = [];
      params.appendTag = '';
      params.mainTag = 'sesspsych';
      params.visible = true;
      params.lineAlpha = 0.5;
      params.angleDivisions = 10;
      params.showIndividualSessions = true;
      params.sessionList = obj.getValidSessions();
      params.verbose = true;
      params = parse_pv_pairs(params, varargin);
      
      if(params.visible)
        visible = 'on';
      else
        visible = 'off';
      end
      
      if(params.verbose)
        slog('Plotting psychometric curves: %s', obj.name);
      end
      
      sessList = params.sessionList(:);
      
      % Use all angles for the list
      validAngles = -diff(abs(obj.data.trialsSummary.ori), [],2); % Left - Right
      angleList = unique(validAngles);
      [~, binnedAngles, angleListIdx] = histcounts(angleList, linspace(-90, 90, params.angleDivisions));
      binnedAnglesCenter = binnedAngles(1:end-1)+diff(binnedAngles)/2;
      hFig = figure('Position', [200 50 1000 300], 'Visible', visible, 'Color', 'w');
      cmap = viridis(length(sessList));
      set(hFig, 'defaultAxesColorOrder', flip(cmap));
      ax = multigap_subplot(1, 5, 'gap_R', 0.01, 'gap_C', [0.05 0.05 0.02 0.02]);
      arrayfun(@(x) hold(x, 'on'), ax);
      arrayfun(@(x) axis(x, 'square'), ax);
      arrayfun(@(x) box(x, 'off'), ax);
      if(params.showIndividualSessions)
        for it0 = sessList
          % Angle difference
          validTrials = ismember(obj.data.trialsSummary.sess, it0);
          choiceDir = obj.data.trialsSummary.choicedir(validTrials);
          outcome = obj.data.trialsSummary.outcome(validTrials);

          performance = nan(size(binnedAnglesCenter));
          performance2 = nan(size(binnedAnglesCenter));
          leftChoice = nan(size(binnedAnglesCenter));
          rightChoice = nan(size(binnedAnglesCenter));
          timeOut = nan(size(binnedAnglesCenter));
          performanceN = nan(size(binnedAnglesCenter));
          performanceN2 = nan(size(binnedAnglesCenter));
          leftChoiceN = nan(size(binnedAnglesCenter));
          rightChoiceN = nan(size(binnedAnglesCenter));
          timeOutN = nan(size(binnedAnglesCenter));
          for it = 1:length(binnedAnglesCenter)
            validAngles = -diff(abs(obj.data.trialsSummary.ori(validTrials, :)), [],2); % Left - Right
            angleSubset = angleList(angleListIdx == it);
            valid = find(arrayfun(@(x)any(x == angleSubset),validAngles));
            performance(it) = sum(outcome(valid) == 1)/length(valid);
            performance2(it) = sum(outcome(valid) == 1)/(sum(outcome(valid) == 1)+sum(outcome(valid) == -1));
            leftChoice(it) = sum(choiceDir(valid) == -1)/length(valid);
            rightChoice(it) = sum(choiceDir(valid) == 1)/length(valid);
            timeOut(it) = sum(choiceDir(valid) == 0)/length(valid);
            performanceN(it) = length(valid);
            performanceN2(it) = (sum(outcome(valid) == 1)+sum(outcome(valid) == -1));
            leftChoiceN(it) = length(valid);
            rightChoiceN(it) = length(valid);
            timeOutN(it) = length(valid);
          end
          cax = ax(1, 1);
          set(hFig, 'currentAxes', cax);
          h = plot(cax, binnedAnglesCenter, performance*100, '-', 'MarkerSize', 4, 'DisplayName', num2str(it0));
          h.Color(4) = params.lineAlpha;
          cax = ax(1, 2);

          set(hFig, 'currentAxes', cax);
          h = plot(cax, binnedAnglesCenter, performance2*100, '-', 'MarkerSize', 4, 'DisplayName', num2str(it0));
          h.Color(4) = params.lineAlpha;

          cax = ax(1, 3);
          set(hFig, 'currentAxes', cax);
          h = plot(cax, binnedAnglesCenter, leftChoice*100, '-', 'MarkerSize', 4, 'DisplayName', num2str(it0));
          h.Color(4) = params.lineAlpha;

          cax = ax(1, 4);
          set(hFig, 'currentAxes', cax);
          h = plot(cax, binnedAnglesCenter, timeOut*100, '-', 'MarkerSize', 4, 'DisplayName', num2str(it0));
          h.Color(4) = params.lineAlpha;

          cax = ax(1, 5);
          set(hFig, 'currentAxes', cax);
          h = plot(cax, binnedAnglesCenter, rightChoice*100, '-', 'MarkerSize', 4, 'DisplayName', num2str(it0));
          h.Color(4) = params.lineAlpha;
        end
      end      
      % Now across all sessoins
      % Angle difference
      validTrials = ismember(obj.data.trialsSummary.sess, sessList);
      choiceDir = obj.data.trialsSummary.choicedir(validTrials);
      outcome = obj.data.trialsSummary.outcome(validTrials);

      performance = nan(size(binnedAnglesCenter));
      performance2 = nan(size(binnedAnglesCenter));
      leftChoice = nan(size(binnedAnglesCenter));
      rightChoice = nan(size(binnedAnglesCenter));
      timeOut = nan(size(binnedAnglesCenter));
      performanceN = nan(size(binnedAnglesCenter));
      performanceN2 = nan(size(binnedAnglesCenter));
      leftChoiceN = nan(size(binnedAnglesCenter));
      rightChoiceN = nan(size(binnedAnglesCenter));
      timeOutN = nan(size(binnedAnglesCenter));
      for it = 1:length(binnedAnglesCenter)
        validAngles = -diff(abs(obj.data.trialsSummary.ori(validTrials, :)), [],2); % Left - Right
        angleSubset = angleList(angleListIdx == it);
        valid = find(arrayfun(@(x)any(x == angleSubset),validAngles));
        performance(it) = sum(outcome(valid) == 1)/length(valid);
        performance2(it) = sum(outcome(valid) == 1)/(sum(outcome(valid) == 1)+sum(outcome(valid) == -1));
        leftChoice(it) = sum(choiceDir(valid) == -1)/length(valid);
        rightChoice(it) = sum(choiceDir(valid) == 1)/length(valid);
        timeOut(it) = sum(choiceDir(valid) == 0)/length(valid);
        performanceN(it) = length(valid);
        performanceN2(it) = (sum(outcome(valid) == 1)+sum(outcome(valid) == -1));
        leftChoiceN(it) = length(valid);
        rightChoiceN(it) = length(valid);
        timeOutN(it) = length(valid);
      end
      
      cax = ax(1, 1);
      set(hFig, 'currentAxes', cax);
      [~, pci] = binofit(round(performance.*performanceN), performanceN);
      sem = diff(pci,[],2)/2;
      m = pci(:,1)+sem;
      h = errorbar(cax, binnedAnglesCenter, m*100,sem*100, '-ok', 'MarkerSize', 4, 'DisplayName', 'global');
      h.MarkerFaceColor =  h.Color;
      
      cax = ax(1, 2);
      set(hFig, 'currentAxes', cax);
      [~, pci] = binofit(round(performance2.*performanceN2), performanceN2);
      sem = diff(pci,[],2)/2;
      m = pci(:,1)+sem;
      h = errorbar(cax, binnedAnglesCenter, m*100,sem*100, '-ok', 'MarkerSize', 4, 'DisplayName', 'global');
      h.MarkerFaceColor =  h.Color;
      
      cax = ax(1, 3);
      set(hFig, 'currentAxes', cax);
      [~, pci] = binofit(round(leftChoice.*leftChoiceN), leftChoiceN);
      sem = diff(pci,[],2)/2;
      m = pci(:,1)+sem;
      h = errorbar(cax, binnedAnglesCenter, m*100,sem*100, '-ok', 'MarkerSize', 4, 'DisplayName', 'global');
      h.MarkerFaceColor =  h.Color;
      
      cax = ax(1, 4);
      set(hFig, 'currentAxes', cax);
      [~, pci] = binofit(round(timeOut.*timeOutN), timeOutN);
      sem = diff(pci,[],2)/2;
      m = pci(:,1)+sem;
      h = errorbar(cax, binnedAnglesCenter, m*100,sem*100, '-ok', 'MarkerSize', 4, 'DisplayName', 'global');
      h.MarkerFaceColor = h.Color;
      
      cax = ax(1, 5);
      set(hFig, 'currentAxes', cax);
      [~, pci] = binofit(round(rightChoice.*rightChoiceN), rightChoiceN);
      sem = diff(pci,[],2)/2;
      m = pci(:,1)+sem;
      h = errorbar(cax, binnedAnglesCenter, m*100,sem*100, '-ok', 'MarkerSize', 4, 'DisplayName', 'global');
      h.MarkerFaceColor = h.Color;

      % Setting axes and labels
      cax = ax(1, 1);
      xlabel(cax, 'Angle difference \langle|\theta_L|-|\theta_R|\rangle');
      ylabel(cax, 'Performance (%)');
      plot(cax, xlim,[50 50],'k:')
      ylim(cax, [0 100]);
      set(cax, 'XTick', linspace(-90,90,5));
      spaceOutAxes(cax);
      
      cax = ax(1, 2);
      xlabel(cax, 'Angle difference \langle|\theta_L|-|\theta_R|\rangle');
      ylabel(cax, 'Performance excl TO (%)');
      plot(cax, xlim,[50 50],'k:')
      ylim(cax, [0 100]);
      set(cax, 'XTick', linspace(-90,90,5));
      spaceOutAxes(cax);
      
      cax = ax(1, 3);
      xlabel(cax, 'Angle difference \langle|\theta_L|-|\theta_R|\rangle');
      title(cax, 'Left');
      ylabel(cax, 'Choice (%)');
      ylim(cax, [0 100]);
      set(cax, 'XTick', linspace(-90,90,5));
      plot(cax, xlim,[50 50],'k:');
      spaceOutAxes(cax);
      
      cax = ax(1, 4);
      xlabel(cax, 'Angle difference \langle|\theta_L|-|\theta_R|\rangle');
      title(cax, 'Time out');
      ylim(cax, [0 100]);
      set(cax, 'YColor', 'none');
      set(cax, 'XTick', linspace(-90,90,5));
      plot(cax, xlim,[50 50],'k:');
      spaceOutAxes(cax);
      
      cax = ax(1, 5);
      xlabel(cax, 'Angle difference \langle|\theta_L|-|\theta_R|\rangle');
      title(cax, 'Right');
      ylim(cax, [0 100]);
      set(cax, 'YAxisLocation','right');
      set(cax, 'XTick', linspace(-90,90,5));
      plot(cax, xlim,[50 50],'k:');
      spaceOutAxes(cax);
      %l = legend(hList);
      %l.Location = 'eastOutside';
      %l.Orientation = 'vertical';
      %l.NumColumns = 2;
      sgtitle(sprintf('Session psychometric curves for: %s', obj.name));

      if(~isempty(params.export))
        obj.exportFig(hFig, params);
      end
      if(~params.visible)
        close(hFig);
      end
    end
    %%
    function [hFig, obj] = plotAverageImage(obj, varargin)
    % PLOTAVERAGEIMAGE Plots the average image
    % USAGE:
    %   obj.plotAverageImage();
    %
    % INPUT arguments:
    %   none
    %
    % INPUT optional arguments ('key' followed by its value):
    %   export    - to automatically export the image to a pdf or png file.
    %               Empty for no export (default). ([], 'pdf', 'png')
    %   mainTag   - main tag to use for the export file name. Default: '_areas' (text)
    %   session   - session data to generate image from (if empty will use
    %   the first valid session) if 0, it will use sessAvgStdImg
    %   appendTag - tag to append to the end of the exported file name. Empty
    %               by default. ([], text)
    %   handleFile - file index with trials to use for the average image (1st one by default)
    %   visible   - if the figure should be visible. Only useful for exports.
    %               If not visible, it will be closed  at the end of the 
    %               script (on by default). ('on','off')
    %   superimposeAreas - none / self /allen / allenLocaNMF
    %   registered
    %   registeredRef
    %
    % OUTPUT arguments:
    %   hFig - Handle of the new figure
    %
    % EXAMPLE:
    %   dataset = preprocessing('animal1', 'animalFile.mat');
    %   dataset = dataset.generateFullData(0, 1, 2.5);
    %   dataset.plotAverageImage('export', 'png');
      
      params.export = [];
      params.mainTag = 'avgImg';
      params.appendTag = '';
      params.visible = true;
      params.handleFile = 1;
      params.superimposeAreas = 'none';
      params.session = [];
      params.verbose = true;
      params.registered = false;
      params.registeredRef = imref2d(size(obj.data.areas));
      params = parse_pv_pairs(params, varargin);

      if(isempty(params.session))
        session = obj.getValidSessions();
        session = session(1);
        data = obj.pullSessionTrials(session);
      else
        if(params.session == 0)
          session = 0;
          data = struct;
          avgNeural = obj.data.sessAvgStdImg;
        else
          session = params.session;
          data = obj.pullSessionTrials(session);
        end
      end
      
      if(params.visible)
        visible = 'on';
      else
        visible = 'off';
      end
      if(params.verbose)
        slog('Plotting session average: %s', obj.name);
      end
      if(session > 0)
        fullNeural = data.DFF;
        %avgNeural = squeeze(nanmean(nanmean(fullNeural, 3), 1));
        %avgNeural = mean(std(fullNeural, [], 3), 2);
        %avgNeural = nanstd(nanmean(fullNeural, 3), [], 2);
        avgNeural = nanmean(nanstd(fullNeural, [], 3), 1);
        %avgNeural = nanstd(nanmean(fullNeural,  3), [], 1);
        %avgNeural = nanstd(nanstd(fullNeural, [], 3), [], 1);
      end      
      [R, C] = size(obj.data.areas);
      avgNeuralImg = reshape(avgNeural, R, C);
      if(params.registered)
        avgNeuralImg = imwarp(avgNeuralImg, obj.data.alignment.transform, 'nearest', 'OutputView', params.registeredRef);
      end
      hFig = figure('Visible', visible, 'Color', 'w');
      %avgNeuralImgNorm = (avgNeuralImg-prctile(avgNeuralImg(:), 15))/(prctile(avgNeuralImg(:), 99)-prctile(avgNeuralImg(:), 15));
      avgNeuralImgNorm = (avgNeuralImg-prctile(avgNeuralImg(:), 1))/(prctile(avgNeuralImg(:), 99)-prctile(avgNeuralImg(:), 1));
      avgNeuralImgNorm(avgNeuralImgNorm > 1) = 1;
      avgNeuralImgNorm(avgNeuralImgNorm < 0) = 0;
      imagesc(avgNeuralImgNorm);
      cmap = viridis(256);
      switch params.superimposeAreas
        case {'self', 'allen', 'allenLocaNMF'}
        hold on;
        if(strcmp(params.superimposeAreas, 'self'))
          areasOrig = obj.data.areas;
          if(params.registered)
            areasOrig = imwarp(areasOrig, obj.data.alignment.transform, 'nearest', 'OutputView', params.registeredRef);
            areasOrig = round(areasOrig);
          end
        elseif(strcmp(params.superimposeAreas, 'allen'))
          areasOrig = obj.data.areasAllen;
        elseif(strcmp(params.superimposeAreas, 'allenLocaNMF'))
          areasOrig = obj.data.areasAllenLocaNMF;
        end
        areasOrigBoundaries = zeros(size(areasOrig));
        for it0 = 1:max(areasOrig(:))
          cArea = zeros(size(areasOrig));
          cArea(areasOrig == it0) = 1;
          CC = bwboundaries(~~cArea,4);
          for it = 1:length(CC)
            for k = 1:size(CC{it}, 1)
              areasOrigBoundaries(CC{it}(k,1),CC{it}(k,2)) = it0;
            end
          end
        end
        h = imagesc(~~areasOrigBoundaries*min(avgNeuralImgNorm(:)));
        h.AlphaData = ~~areasOrigBoundaries;
        box on;
        colorbar;
        %caxis([0.5 length(cmap)+0.5]);
        hold on;
        for it = 1:max(areasOrig(:))
          [r, c] = find(areasOrig == it);
          text(mean(c)-1, mean(r), num2str(it),'Color', 'k', 'FontWeight','bold');
        end
        otherwise
      end
      if(params.registered)
        title(sprintf('Avg Img for: %s session: %d - registered', obj.name, session));
      else
        title(sprintf('Avg Img for: %s session: %d', obj.name, session));
      end
      colormap(cmap);
      axis equal tight
      
      if(~isempty(params.export))
        obj.exportFig(hFig, params);
      end
      if(~params.visible)
        close(hFig);
      end
      obj.data.avgImg = avgNeuralImgNorm;
    end
    %%
    function hFig = plotSessionChoice(obj, varargin)
    % PLOTSESSIONCHOICE Plots the choices for each session
    % USAGE:
    %   obj.plotSessionChoice();
    %
    % INPUT arguments:
    %   none
    %
    % INPUT optional arguments ('key' followed by its value):
    %  doBiasPlot - Plot the bias instead (true / false (default))
    %   export    - to automatically export the image to a pdf or png file.
    %               Empty for no export. ([] (default), 'pdf', 'png')
    %   mainTag   - main tag to use for the export file name. (txt, '_fps' (default))
    %   appendTag - tag to append to the end of the exported file name. Empty
    %               by default. ([] (default), text)
    %   visible   - if the figure should be visible. Only useful for exports.
    %               If not visible, it will be closed  at the end of the 
    %               script (on by default). (true/false)
    %   verbose   - displays additional info (true/false)
    %
    % OUTPUT arguments:
    %   hFig - Handle of the new figure
    %
    % EXAMPLE:
    %   dset = wf.dataset.load(pwd, 'animal1');
    %   dset.plotSessionChoice('export', 'pdf');
      
      params.doBiasPlot = false;
      params.export = [];
      params.appendTag = '';
      params.mainTag = 'sesschoice';
      params.verbose = true;
      params.visible = true;
      params = parse_pv_pairs(params, varargin);
      
      if(params.visible)
        visible = 'on';
      else
        visible = 'off';
      end
      
      if(params.verbose)
        slog('Plotting session choices: %s', obj.name);
      end
      
      sessList = obj.getValidSessions();
      
      sessLRTO = zeros(max(obj.data.trialsSummary.sess), 3);
      for it0 = sessList
        sessLRTO(it0, :) = [sum(obj.data.trialsSummary.choicedir(ismember(obj.data.trialsSummary.sess, it0)) == -1), ...
                            sum(obj.data.trialsSummary.choicedir(ismember(obj.data.trialsSummary.sess, it0)) ==  1), ...
                            sum(obj.data.trialsSummary.choicedir(ismember(obj.data.trialsSummary.sess, it0)) ==  0)]/ sum(ismember(obj.data.trialsSummary.sess, it0));
      end
      
      hFig = figure('Visible', visible, 'Color', 'w');
      hold on;
      if(params.doBiasPlot)
        bar(sessLRTO(:,1)-sessLRTO(:,2), 0.8);
%         bb = bar(sessLRTO(:,3)-sessLRTO(:,3)/2, 0.4);
%         bb.FaceAlpha = 0.5;
%         bb.EdgeAlpha = 0.5;
%         bb2 = bar(-(sessLRTO(:,3)-sessLRTO(:,3)/2), 0.4);
%         bb2.FaceColor = bb.FaceColor;
%         bb2.EdgeColor = bb.EdgeColor;
%         bb2.FaceAlpha = 0.5;
%         bb2.EdgeAlpha = 0.5;
        maxD = max(abs(sessLRTO(:,1)-sessLRTO(:,2)));
        ylim([-1.1 1.1]*maxD);
        title(sprintf('Session choice bias: %s', obj.name));
        ylabel('choice bias fraction (Left - Right)');
      else
        bar(sessLRTO, 0.8, 'stacked');
        l = legend('L', 'R', 'TO');
        ylim([0 1]);
        title(sprintf('Session choice: %s', obj.name));
        ylabel('choice fraction');
      end
      xlim([1 max(sessList)]+[-0.4 0.4]);

      xlabel('Session');
      box on;

      l.Orientation = 'horizontal';
      l.Location = 'northOutside';
      
      l.Box = 'off';
      
      if(~isempty(params.export))
        obj.exportFig(hFig, params);
      end
      if(strcmp(params.visible, 'off'))
        close(hFig);
      end
    end
    %%
    function hFig = plotAreas(obj, varargin)
    % PLOTAREAS Plots the retinotopic areas
    % USAGE:
    %   obj.plotAreas();
    %
    % INPUT arguments:
    %   none
    %
    % INPUT optional arguments ('key' followed by its value):
    %   export    - to automatically export the image to a pdf or png file.
    %               Empty for no export (default). ([], 'pdf', 'png')
    %   mainTag   - main tag to use for the export file name. Default: '_areas' (text)
    %   appendTag - tag to append to the end of the exported file name. Empty
    %               by default. ([], text)
    %   visible   - if the figure should be visible. Only useful for exports.
    %               If not visible, it will be closed  at the end of the 
    %               script (on by default). ('on','off')
    %   onlyContours - will plot only the contours of the areas (false by default). (true, false)
    %   handle - axes handle to plot the figure at. If empty, will create a new figure (default)
    %   type - original / registered / allen / loca
    %   registeredRef - imref2d object
    %   verbose   - displays additional info (true/false)
    % 
    %
    % OUTPUT arguments:
    %   hFig - Handle of the new figure
    %
    % EXAMPLE:
    %   dset = wf.dataset.load(pwd, 'animal1');
    %   dset.plotAreas('export', 'pdf');
      
      params.rounding = [];
      params.export = [];
      params.mainTag = 'areas';
      params.appendTag = '';
      params.visible = true;
      params.onlyContours = false;
      params.handle = [];
      params.type = 'original';
      params.registeredRef = imref2d(size(obj.data.areas));
      params.verbose = true;
      params = parse_pv_pairs(params, varargin);
      
      areasOrig = obj.data.areas;
      switch params.type
        case 'original'
        case 'registered'
          areasOrig = imwarp(areasOrig, obj.data.alignment.transform, 'nearest', 'OutputView', params.registeredRef);
          areasOrig = round(areasOrig);
        case 'allen'
          %params.registeredRef = imref2d(size(obj.data.areasAllen));
          %areasOrig = imwarp(obj.data.areasAllen, obj.data.alignment.transform, 'nearest', 'OutputView', params.registeredRef);
          areasOrig = obj.data.areasAllen;
        case 'loca'
          %params.registeredRef = imref2d(size(obj.data.areasAllen));
          %areasOrig = imwarp(obj.data.areasAllenLocaNMF, obj.data.alignment.transform, 'nearest', 'OutputView', params.registeredRef);
          areasOrig = obj.data.areasAllenLocaNMF;
      end
      cmap = viridis(max(areasOrig(:)));
      if(isempty(params.handle))
        hFig = figure('Visible', params.visible, 'Color', 'w');
      else
        hFig = params.handle;
      end
      if(params.onlyContours)
        areasOrigBoundaries = zeros(size(areasOrig));
        for it0 = 1:max(areasOrig(:))
          cArea = zeros(size(areasOrig));
          cArea(areasOrig == it0) = 1;
          CC = bwboundaries(~~cArea,4);
          for it = 1:length(CC)
            for k = 1:size(CC{it}, 1)
              areasOrigBoundaries(CC{it}(k,1),CC{it}(k,2)) = it0;
            end
          end
        end
        h = imagesc(areasOrigBoundaries);
        h.AlphaData = ~~areasOrigBoundaries;
      else
        areasOrigBoundaries = zeros(size(areasOrig));
        CC = bwboundaries(~~areasOrig,4);
        for it = 1:length(CC)
          for k = 1:size(CC{it}, 1)
            areasOrigBoundaries(CC{it}(k,1),CC{it}(k,2)) = it;
          end
        end
        h = imagesc(areasOrig);
        h.AlphaData = ~~areasOrig;
        hold on;
        h = imagesc(~~areasOrigBoundaries);
        h.AlphaData = ~~areasOrigBoundaries*0.15;
      end
      
      colormap(cmap);
      axis equal tight;
      switch params.type
        case 'registered'
          title(sprintf('Areas for: %s - registered', obj.name));
        case 'allen'
          title(sprintf('Areas for: %s - allen', obj.name));
        case 'loca'
          title(sprintf('Areas for: %s - loca', obj.name));
        case 'original'
          title(sprintf('Areas for: %s', obj.name));
      end
      box on;
      cb = colorbar;
      caxis([0.5 length(cmap)+0.5]);
      hold on;
      try
        if(params.registered)
          hh = plot(obj.data.alignment.triangle, 'EdgeColor', [1 1 1], 'FaceColor', [1 1 1], 'LineWidth', 1);
        else
          hh = plot(obj.data.retinotopyCoordinates.alignmentTriangle, 'EdgeColor', [1 1 1], 'FaceColor', [1 1 1], 'LineWidth', 1);
        end
        hh.FaceAlpha = 0.5;
        hh.EdgeAlpha = 0.5;
      end
      areaNames = string(1:max(areasOrig(:)));
      for it = 1:max(areasOrig(:))
        validID = false;
        [r, c] = find(areasOrig == it);
        %dset.data.retinotopyCoordinates.areaID
        try
          fnames = fieldnames(obj.data.retinotopyCoordinates.areaID);
          for it2 = 1:length(fnames)
            if(obj.data.retinotopyCoordinates.areaID.(fnames{it2}) == it)
              validID = true;
              text(mean(c)-1, mean(r), fnames{it2}, 'Color', 'k', 'FontWeight','bold');
              areaNames{it} = fnames{it2};
            end
          end
        end
        if(~validID)
          text(mean(c)-1, mean(r), num2str(it),'Color', 'k', 'FontWeight','bold');
        end
        %text(mean(c)-1, mean(r), num2str(it),'Color', 'k');
      end
      cb.Ticks = 1:max(areasOrig(:));
      cb.TickLabels = areaNames;
      switch params.type
        case 'original'
        otherwise
        xlim(params.registeredRef.XWorldLimits);
        ylim(params.registeredRef.YWorldLimits);
      end
      if(~isempty(params.export))
        obj.exportFig(hFig, params);
      end
      if(strcmp(params.visible, 'off'))
        close(hFig);
      end
    end
    
    %% --------------------------------------------------------------------
    %%% Preprocessing methods
    %%% -------------------------------------------------------------------
    
    function obj = generateTimeSeries(obj, preStimulusStart, openLoopStart, closedLoopStart, closedLoopEnd, varargin)
    % GENERATETIMESERIES Creates a single file will all info from the trials
    % USAGE:
    %   obj.generateFullData();
    %
    % INPUT arguments:
    %   preStimulusStart - time (in sec) when the prestimulus starts in a trial
    %   openLoopStart    - time (in sec) when the open loop starts in a trial
    %   closedLoopStart  - time (in sec) when the closed loop starts in a trial
    %   closedLoopEnd    - time (in sec) when the closed loop ends in a trial
    %
    % INPUT optional arguments ('key' followed by its value):
    %   mainTag   - main tag to use for the export file name. Default: '_fullData' (text)
    %   appendTag - tag to append to the end of the exported file name. Empty
    %               by default. ([], text)
    %   bar       - true to show the progressbar (true (default) / false)
    %   fps       - frame rate to resample the data to. If empty, will use the one found on the first trial
    %   lowPassFilter - if a low pass filter should be used (zero-phase) on the data. Empty for no filter. Or frequency in Hz
    %   frameRateRoundingDecimals - will round the framerate value to the nearest Xth decimal (default X=1, empty for no rounding)
    %
    % OUTPUT arguments:
    %   obj - class object
    %
    % EXAMPLE:
    %   dset = wf.dataset.load(pwd, 'animal1');
    %   dset = dset.generateTimeSeries(0, 1, 2.5, 12.5);
       
      params.mainTag = 'timeSeries';
      params.appendTag = '';
      params.lowPassFilter = 8;
      params.bar = true;
      params.fps = [];
      params.restBaseline = false;
      params.pixelBaseline = false;
      params.register = false;
      params.frameRateRoundingDecimals = 1;
      params = parse_pv_pairs(params, varargin);
      
      sessList = obj.getValidSessions();
      validTrials = find(ismember(obj.data.trialsSummary.sess, sessList));
      totalTrials = length(obj.data.trialsSummary.sess);
      
      fullData = struct;

      fullData.validTrials = zeros(totalTrials, 1);
      fullData.validTrials(validTrials) = 1;
      fullData.handleFile = zeros(totalTrials, 1);
      fullData.invalidTrials = ones(totalTrials, 1);
      fullData.invalidTrials(validTrials) = 0;
      % Now remove trials with sMo
      if(isfield(obj.data.trialsSummary, 'sMo'))
        fullData.validTrials(obj.data.trialsSummary.sMo(:,1) ~= 0) = 0;
        fullData.invalidTrials(obj.data.trialsSummary.sMo(:,1) ~= 0) = 1;
        fprintf(sprintf('Removed %d trials from stimulus motion\n', sum(~~obj.data.trialsSummary.sMo(:,1))));
      end
      validTrials = find(fullData.validTrials);
      invalidTrialCount = 0;
      fullData.params = struct;
      fullData.params.openLoopStart = openLoopStart;
      fullData.params.closedLoopStart = closedLoopStart;
      fullData.params.preStimulusStart = preStimulusStart;
      fullData.params.closedLoopEnd = closedLoopEnd;
      fullData.params.lowPassFilter = params.lowPassFilter;
      fullData.params.restBaseline = params.restBaseline;
      fullData.params.pixelBaseline = params.pixelBaseline;
      fullData.params.selAnimal = obj.name;
      fullData.params.sessList = sessList;
      fullData.params.fps = params.fps;
      
      fullData.decisionTime = nan(totalTrials, 1);
      fullData.decisionTimeCorr = nan(totalTrials, 1);
      fullData.lowPassFilter = [];
      
      if(params.bar)
        ncbar(sprintf('Preprocessing trials: %s', obj.name));
      end

      validTrialsPerSession = cell(length(sessList), 1);
      
      for it0 = 1:length(sessList)
        validTrialsPerSession{it0} = validTrials(find(obj.data.trialsSummary.sess(validTrials) == sessList(it0)));
      end

      trialOffset = 0;
      handlesList = [];
         
      for it0 = 1:length(validTrials)
        curTrial = validTrials(it0);
        curSession = obj.data.trialsSummary.sess(curTrial);
        curSessionIdx = find(sessList == curSession);
        
        W = obj.pullSingleTrial(curTrial);
        try
          wheelPos = W.wheelpos;
          pupilPosX = W.ppos(1,:);
          pupilPosY = W.ppos(2,:);
          pupilArea = W.parea;

          if(~isempty(params.frameRateRoundingDecimals))
            fps = round(W.fps*10^params.frameRateRoundingDecimals)/10^params.frameRateRoundingDecimals;
          end
          timeNeural = ((1:size(W.resp,3))-1)/fps(1);
          timeBehavior = ((1:size(W.parea,2))-1)/fps(2);

          % Set the correct framerate
          if(it0 == 1)
            if(isempty(params.fps))
              fullData.fps = fps(1);
            else
              fullData.fps = params.fps;
            end
            % Set all the time vectors
            fpsTarget = fullData.fps;
            timeTarget = preStimulusStart:1/fpsTarget:closedLoopEnd;
            firstPrestimulus = find(timeTarget >= preStimulusStart, 1, 'first');
            lastPrestimulus = find(timeTarget < openLoopStart, 1, 'last');
            firstOpenLoop = find(timeTarget >= openLoopStart, 1, 'first');
            lastOpenLoop = find(timeTarget < closedLoopStart, 1, 'last');
            firstClosedLoop = find(timeTarget >= closedLoopStart, 1, 'first');
            lastClosedLoop = find(timeTarget < closedLoopEnd, 1, 'last');

            timePrestimulus = timeTarget(firstPrestimulus:lastPrestimulus);
            timeOpenLoop = timeTarget(firstOpenLoop:lastOpenLoop);
            timeClosedLoop = timeTarget(firstClosedLoop:lastClosedLoop);
            timeTarget = timeTarget(firstPrestimulus:lastClosedLoop);
            
            fullData.timePrestimulus = timePrestimulus;
            fullData.timeOpenLoop = timeOpenLoop;
            fullData.timeClosedLoop = timeClosedLoop; % Not to use
            fullData.firstPrestimulus = firstPrestimulus;
            fullData.lastPrestimulus = lastPrestimulus;
            fullData.firstOpenLoop = firstOpenLoop;
            fullData.lastOpenLoop = lastOpenLoop;
            fullData.firstClosedLoop = firstClosedLoop;
            fullData.lastClosedLoop = ceil(closedLoopEnd*fullData.fps);
            fullData.timeSeries = cell(3, 1);
            for it1 = 1:3 % pre/OL/CL
              switch it1
                case 1
                  selFrames = fullData.firstPrestimulus:fullData.lastPrestimulus;
                case 2
                  selFrames = fullData.firstOpenLoop:fullData.lastOpenLoop;
                case 3
                  
                  selFrames = fullData.firstClosedLoop:fullData.lastClosedLoop;
              end
              % Changing to column-major order. So data from each trial are together in the file
              fullData.timeSeries{it1}.pX = nan(length(selFrames), length(validTrialsPerSession{curSessionIdx}));
              fullData.timeSeries{it1}.pY = nan(length(selFrames), length(validTrialsPerSession{curSessionIdx}));
              fullData.timeSeries{it1}.pA = nan(length(selFrames), length(validTrialsPerSession{curSessionIdx}));
              fullData.timeSeries{it1}.pS = nan(length(selFrames), length(validTrialsPerSession{curSessionIdx}));
              fullData.timeSeries{it1}.wP = nan(length(selFrames), length(validTrialsPerSession{curSessionIdx}));
              fullData.timeSeries{it1}.wV = nan(length(selFrames), length(validTrialsPerSession{curSessionIdx}));
              fullData.timeSeries{it1}.T  = nan(length(selFrames), length(validTrialsPerSession{curSessionIdx}));
              fullData.timeSeries{it1}.indexRaw  = nan(1, length(validTrialsPerSession{curSessionIdx}));
              fullData.timeSeries{it1}.indexValid = nan(1, length(validTrialsPerSession{curSessionIdx}));
              fullData.timeSeries{it1}.neural = nan(size(W.resp, 2), length(selFrames), length(validTrialsPerSession{curSessionIdx}));
            end
            % Low pass filter check
            if(~isempty(params.lowPassFilter))
              lpFilt = designfilt('lowpassiir', 'FilterOrder', 6, ...
               'PassbandFrequency', params.lowPassFilter,'PassbandRipple', 0.1, ...
               'SampleRate', fullData.fps);
            else
              lpFilt = [];
            end
            fullData.lowPassFilter = lpFilt;
          end
          % Interpolate data as needed
          if(fps(1) ~= fpsTarget)
              fprintf('Warning, diff frame rate on trial %d (%.2f). Will interpolate to %.2f\n', curTrial, fps(1), fpsTarget);
              interpMovie = true;
          end
          % X pupil
          %pupilPosX = interp1(timeBehavior, pupilPosX, timeTarget, 'pchip', 'extrap');
          pupilPosX = interp1(timeBehavior, pupilPosX, timeTarget, 'pchip', NaN);
          % Y pupil
          pupilPosY = interp1(timeBehavior, pupilPosY, timeTarget, 'pchip', NaN);
          % Area pupil
          pupilPosA = interp1(timeBehavior, pupilArea, timeTarget, 'pchip', NaN);
          % Saccade signature
          saccadesT = zeros(size(timeTarget));
          saccadeTimes = timeBehavior(W.pmvton);
          for its = 1:length(W.pmvton)
            [~, valid] = min(abs(saccadeTimes(its)-timeTarget));
            saccadesT(valid) = 1;
          end
          % Wheel position 
          wheelP = interp1(wheelPos(:,2), wheelPos(:,1), timeTarget, 'pchip', NaN);
          % Wheel position Diff
          wheelV = interp1(wheelPos(1:end-1,2), diff(wheelPos(:,1)), timeTarget, 'pchip', NaN);
          % Wheel fix - These sessions had a different gain on the wheel encoder
          switch obj.name
            case 'A15100'
              if(obj.data.trialsSummary.sess(curTrial) <= 9)
                wheelP = wheelP/6;
                wheelV = wheelV/6;
              end
            case 'A15098'
              if(obj.data.trialsSummary.sess(curTrial) <= 18)
                wheelP = wheelP/6;
                wheelV = wheelV/6;
              end
          end
          movDataOrig = squeeze(W.resp(1, :, :))';

          if(any(isnan(movDataOrig(:))))
            invalidFrames = find(isnan(sum(movDataOrig, 2)));
            fprintf('Warning. Found NaN on trial: %d - Frames:', curTrial);
            fprintf(' %d', invalidFrames);
            fprintf('\n');
          end
          %if(interpMovie)
          validFrames = find(~isnan(sum(movDataOrig, 2)));
          validFrames(validFrames > length(timeTarget)) = [];
          if(all(abs(timeNeural(validFrames)-timeTarget(validFrames)) < 1e-5))
            movData = nan(length(timeTarget), size(movDataOrig, 2));
            movData(validFrames, :) = movDataOrig(validFrames, :);
          else
            movData = interp1(timeNeural(validFrames), movDataOrig(validFrames, :), timeTarget, 'pchip', NaN);
          end
          %end

          try
            fullData.decisionTime(curTrial) = W.eventtimes.responseMadeTime-W.eventtimes.interactiveStartedTime;
          catch
            fullData.decisionTime(curTrial) = nan;
          end
          fullData.decisionTimeCorr(curTrial) = timeNeural(end);
          
          % Low pass filtering (exclude nans)
          if(~isempty(lpFilt))
            validRange = find(~isnan(mean(movData,2)));
            movData(validRange, :) = filtfilt(lpFilt, movData(validRange, :));
          end
          
          baseF = [];
          % 0-centering on the resting stage
          if(params.restBaseline)
            baseF = mean(mean(movData(firstPrestimulus:lastPrestimulus, :)));
            movData = movData - baseF;
          end
          % Using the first 3 frames
          if(params.pixelBaseline)
            baseF = mean(movData(firstPrestimulus:(firstPrestimulus+2), :), 1)';
            movData = movData - repmat(baseF, [1 size(movData, 1)])';
          end
          % Split the frames
          for it1 = 1:3
            switch it1
              case 1
                selFrames = firstPrestimulus:lastPrestimulus;
              case 2
                selFrames = firstOpenLoop:lastOpenLoop;
              case 3
                selFrames = firstClosedLoop:lastClosedLoop;
            end
            fullData.timeSeries{it1}.pX(1:length(selFrames), it0-trialOffset) = pupilPosX(selFrames);
            fullData.timeSeries{it1}.pY(1:length(selFrames), it0-trialOffset) = pupilPosY(selFrames);
            fullData.timeSeries{it1}.pA(1:length(selFrames), it0-trialOffset) = pupilPosA(selFrames);
            fullData.timeSeries{it1}.pS(1:length(selFrames), it0-trialOffset) = saccadesT(selFrames)';
            fullData.timeSeries{it1}.wP(1:length(selFrames), it0-trialOffset) = wheelP(selFrames);
            fullData.timeSeries{it1}.wV(1:length(selFrames), it0-trialOffset) = wheelV(selFrames);
            fullData.timeSeries{it1}.T(1:length(selFrames), it0-trialOffset) = timeTarget(selFrames);
            fullData.timeSeries{it1}.indexRaw(it0-trialOffset) = curTrial;
            fullData.timeSeries{it1}.indexValid(it0-trialOffset) = it0;
            fullData.timeSeries{it1}.neural(:, 1:length(selFrames), it0-trialOffset) = movData(selFrames, :)';
          end
        catch ME
          invalidTrialCount = invalidTrialCount + 1;
          fullData.invalidTrials(curTrial) = 1;
          fullData.validTrials(curTrial) = 0;
          fprintf('Invalid trial: %d\n', curTrial)
          fprintf('%s\n', ME.message);
        end
        if(validTrialsPerSession{curSessionIdx}(end) == curTrial)
          %curSavePoint = find(savePoints == it0);
          %savePointsRangeText = sprintf('_trials_%dto%d', savePoints(curSavePoint-1)+1, savePoints(curSavePoint));
          savePointsRangeText = sprintf('trials_sess_%d', sessList(curSessionIdx));
          fprintf('Saving time series from trials from session %d (%d/%d)\n', sessList(curSessionIdx), it0, length(validTrials));
          saveFile1 = sprintf('%s_%s_preStimulus_%s%s.mat', obj.name, params.mainTag, savePointsRangeText, params.appendTag);
          fullDataFile = fullfile(obj.rootFolder, obj.name, obj.dataFolder, saveFile1);
          timeSeries = fullData.timeSeries{1};
          save(fullDataFile, '-struct', 'timeSeries');

          saveFile2 = sprintf('%s_%s_openLoop_%s%s.mat', obj.name, params.mainTag, savePointsRangeText, params.appendTag);
          fullDataFile = fullfile(obj.rootFolder, obj.name, obj.dataFolder, saveFile2);
          timeSeries = fullData.timeSeries{2};
          save(fullDataFile, '-struct', 'timeSeries');

          saveFile3 = sprintf('%s_%s_closedLoop_%s%s.mat', obj.name, params.mainTag, savePointsRangeText, params.appendTag);
          fullDataFile = fullfile(obj.rootFolder, obj.name, obj.dataFolder, saveFile3);
          timeSeries = fullData.timeSeries{3};
          save(fullDataFile, '-struct', 'timeSeries');
          
          fullData.handleFile(validTrialsPerSession{curSessionIdx}) = curSessionIdx;
          handlesList = [handlesList; {saveFile1, saveFile2, saveFile3}];
          if(it0 < length(validTrials))
            trialOffset = it0;
            for it1 = 1:3
              switch it1
                case 1
                  selFrames = firstPrestimulus:lastPrestimulus;
                case 2
                  selFrames = firstOpenLoop:lastOpenLoop;
                case 3
                  selFrames = firstClosedLoop:lastClosedLoop;
              end
              
              fullData.timeSeries{it1}.pX = nan(length(selFrames), length(validTrialsPerSession{curSessionIdx+1}));
              fullData.timeSeries{it1}.pY = nan(length(selFrames), length(validTrialsPerSession{curSessionIdx+1}));
              fullData.timeSeries{it1}.pA = nan(length(selFrames), length(validTrialsPerSession{curSessionIdx+1}));
              fullData.timeSeries{it1}.pS = nan(length(selFrames), length(validTrialsPerSession{curSessionIdx+1}));
              fullData.timeSeries{it1}.wP = nan(length(selFrames), length(validTrialsPerSession{curSessionIdx+1}));
              fullData.timeSeries{it1}.wV = nan(length(selFrames), length(validTrialsPerSession{curSessionIdx+1}));
              fullData.timeSeries{it1}.T  = nan(length(selFrames), length(validTrialsPerSession{curSessionIdx+1}));
              fullData.timeSeries{it1}.indexRaw  = nan(1, length(validTrialsPerSession{curSessionIdx+1}));
              fullData.timeSeries{it1}.indexValid  = nan(1, length(validTrialsPerSession{curSessionIdx+1}));
              fullData.timeSeries{it1}.neural = nan(size(W.resp, 2), length(selFrames), length(validTrialsPerSession{curSessionIdx+1}));
            end
          end
        end
        if(params.bar)
          ncbar.update(it0/length(validTrials));
        end
      end
      if(params.bar)
        ncbar.close();
      end
      if(invalidTrialCount > 0)
        fprintf('Found %d invalid trials\n', invalidTrialCount);
      end
      fullData.validTrials = find(fullData.validTrials);
      fullData.invalidTrials = find(fullData.invalidTrials);
      fullData.baseF = baseF;
      %%%
      saveFile = sprintf('%s_%s_Info%s.mat', obj.name, params.mainTag, params.appendTag);
      fullDataFile = fullfile(obj.rootFolder, obj.name, obj.dataFolder, saveFile);
      fullData = rmfield(fullData, 'timeSeries');
      fullData.handles = handlesList;
      obj.data.timeSeriesInfo = fullData;
      obj.data.timeSeriesInfo.params = params;
      timeSeriesInfo = fullData; %#ok<PROPLC>
      save(fullDataFile, 'timeSeriesInfo');
      fprintf('Saving done!\n');
      obj.data.timeSeriesInfoHandle = fullDataFile;
    end

    %%
    function data = pullSingleTrial(obj, trialNumber)
    % PULLSINGLETRIAL Pulls the data from a single trial (the W file)
    % USAGE:
    %   data = obj.pullSingleTrial(trialNumber);
    %
    % INPUT arguments:
    %   trialNumber - index of the trial
    %
    % OUTPUT arguments:
    %   data - data from the trial
    %
    % EXAMPLE:
    %   dset = wf.dataset.load(pwd, 'animal1');
    %   data = dset.pullSingleTrial(21);
      try
        data = load(fullfile(obj.info.trialsSummaryFolder, 'tile', [num2str(trialNumber) '.mat']));
      catch
        data = [];
        slog.warning('Could not find trial: %d', num2str(trialNumber));
      end
    end
    %%
    function data = pullSessionTrials(obj, sessionNumber, varargin)
    % PULLSESSIONTRIALS Pulls all the trials from a given session (concatenated)
    % USAGE:
    %   data = obj.pullSessionTrials(sessionNumber);
    %
    % INPUT arguments:
    %   sessionNumber - index of the session
    %
    % INPUT optional arguments ('key' followed by its value):
    %   force   - to force generating the dataset / rather than loading
    %   preexisting stuff
    %   fps     - if not empty, will rescale all data to that fps
    %   lowPassFilter  - if not empty, will apply a lowpass filter at the desire frequency to smooth the
    %   type       - empty / normal / raw / clean
    %   neural data
    %  frameRateRoundingDecimals
    %  maxTime
    %
    % OUTPUT arguments:
    %   data - data from the trials
    %
    % EXAMPLE:
    %   dset = wf.dataset.load(pwd, 'animal1');
    %   data = dset.pullSessionTrials(21);
      
      params.bar = false;
      params.force = false;
      params.fps = [];
      params.frameRateRoundingDecimals = 1;
      params.lowPassFilter = 8;
      params.type = [];
      params.maxTime = [];
      
      params = parse_pv_pairs(params, varargin);
      %fullfile(obj.rootFolder, obj.name, obj.dataFolder, sprintf('%s_session_%d.mat', obj.name, sessionNumber))
      %exist(fullfile(obj.rootFolder, obj.name, obj.dataFolder, sprintf('%s_session_%d.mat', obj.name, sessionNumber)), 'file')
      if(~params.force && exist(fullfile(obj.rootFolder, obj.name, obj.dataFolder, sprintf('%s_session_%d.mat', obj.name, sessionNumber)), 'file'))
        if(isempty(params.type))
          data = load(fullfile(obj.rootFolder, obj.name, obj.dataFolder, sprintf('%s_session_%d.mat', obj.name, sessionNumber)));
        else
          switch params.type
            case 'normal'
              data = load(fullfile(obj.rootFolder, obj.name, obj.dataFolder, sprintf('%s_session_%d.mat', obj.name, sessionNumber)));
            case 'raw'
              data = load(fullfile(obj.rootFolder, obj.name, obj.dataFolder, sprintf('%s_sessionRaw_%d.mat', obj.name, sessionNumber)));
            case 'clean'
              data = load(fullfile(obj.rootFolder, obj.name, obj.dataFolder, sprintf('%s_sessionClean_%d.mat', obj.name, sessionNumber)));
          end
        end
        slog('Session %d data succesfully loaded from file: %s', sessionNumber, fullfile(obj.rootFolder, obj.name, obj.dataFolder, sprintf('%s_session_%d.mat', obj.name, sessionNumber)));
        return;
      end
      if(params.bar)
        ncbar(sprintf('Preprocessing trials for session: %d', sessionNumber));
      end
      validTrials = find(sum(obj.data.trialsSummary.sess == sessionNumber, 2));
      data = struct;
      fnames = {'parea', 'pax', 'pfix', 'pmvton', 'ppos', 'resp', 'wheelmvtacc', 'wheelmvton', 'wheelmvtsz', 'wheelpos'};
      for it2 = 1:length(fnames)
        data.(fnames{it2}) = [];
      end
      data.timeF = []; % Extra
      data.timeP = []; % Extra
      data.timeReal = []; % Extra
      data.trialIDF = []; % Extra
      data.trialIDP = []; % Extra
      data.trialIDW = []; % Extra
      data.trialIDWM = []; % Extra
      data.trialIDS = []; % Extra
      

      interpMovie = false;
      for it1 = 1:length(validTrials)
        curTrial = validTrials(it1);
        trialData = obj.pullSingleTrial(curTrial);
        if(~isempty(params.frameRateRoundingDecimals))
          fps = round(trialData.fps*10^params.frameRateRoundingDecimals)/10^params.frameRateRoundingDecimals;
        end
        if(it1 == 1) % Need to load a trial to get this data first
          data.fps = params.fps;
          if(isempty(data.fps))
            data.fps = fps(1); % The neural data
          end
          if(~isempty(params.lowPassFilter))
            lpFilt = designfilt('lowpassiir', 'FilterOrder', 6, ...
               'PassbandFrequency', params.lowPassFilter,'PassbandRipple', 0.1, ...
               'SampleRate', data.fps);
          else
            lpFilt = [];
          end
          data.lpFilt = lpFilt;
          if(fps(1) ~= data.fps)
              fprintf('Warning, diff frame rate on F data (%.2f). Will interpolate to %.2f\n', fps(1), data.fps);
              interpMovie = true;
          end
        end
        
        for it2 = 1:length(fnames)
          switch fnames{it2}
            case 'pmvton'
              data.(fnames{it2}) = [data.(fnames{it2}), trialData.(fnames{it2})'/fps(2)]; % It's in pupil index coordinates
            case {'wheelmvtacc', 'wheelmvton', 'wheelmvtsz', 'wheelpos'}
              %trialData.(fnames{it2})
              data.(fnames{it2}) = [data.(fnames{it2}), trialData.(fnames{it2})'];
            case 'resp'
              nData = squeeze(trialData.(fnames{it2}));
              data.(fnames{it2}) = [data.(fnames{it2}), nData];
              data.timeF = [data.timeF, ((1:size(nData, 2))-1)/fps(1)];
              try
                fID = fopen(trialData.resppix.fnm,'r');
                fread(fID,3, '*uint32'); % image data size
                posixTime = fread(fID,1,'*double'); % UNIX time stamp
                fclose(fID);
              catch
                posixTime = 0;
              end
              %data.timeReal = [data.timeReal, datetime(posixTime+((1:size(nData, 2))-1)/trialData.fps(1), 'ConvertFrom', 'posixtime')];
              data.timeReal = [data.timeReal, posixTime+((1:size(nData, 2))-1)/fps(1)];
            otherwise
              data.(fnames{it2}) = [data.(fnames{it2}), trialData.(fnames{it2})];
          end
        end
        data.timeP = [data.timeP, ((1:size(trialData.parea, 2))-1)/fps(2)];
        data.trialIDF = [data.trialIDF, curTrial*ones(1, size(nData, 2))];
        data.trialIDP = [data.trialIDP, curTrial*ones(1, size(trialData.parea, 2))];
        data.trialIDW = [data.trialIDW, curTrial*ones(1, size(trialData.wheelpos, 1))];
        data.trialIDWM = [data.trialIDWM, curTrial*ones(1, size(trialData.wheelmvtacc, 1))];
        data.trialIDS = [data.trialIDS, curTrial*ones(1, size(trialData.pmvton, 1))];
        if(params.bar)
          ncbar.update(it1/length(validTrials));
        end
      end

      % Now it's time to reshape the data intro trials
      % Nan'ing
      reshapedData = struct;
      trialID = unique(data.trialIDP);
      Ntrials = length(trialID);
      reshapedData.trialID = trialID;
      % Start with the pupil
      maxFrames = max(hist(data.trialIDP, 1:max(trialID)));
      timePraw = ((1:maxFrames)-1)/fps(2);
      reshapedData.timeP = 0:1/data.fps:max(timePraw);
      maxFrames = length(reshapedData.timeP);
      
      fnames = {'pupilArea', 'pupilX', 'pupilY', 'saccades'};
      for it1 = 1:length(fnames)
        reshapedData.(fnames{it1}) = nan(maxFrames, Ntrials);
      end
      
      % Do Fluorescence
      maxFrames = max(hist(data.trialIDF, 1:max(trialID)));
      timeFraw = ((1:maxFrames)-1)/fps(1);
      reshapedData.timeF = 0:1/data.fps:max(timeFraw);
      maxFrames = length(reshapedData.timeF);
      reshapedData.DFF = nan(maxFrames, size(data.resp, 1), Ntrials);

      % Do the wheel
      [maxFrames, maxFramesIDW] = max(hist(data.trialIDW, 1:max(trialID)));
      timeWraw = data.wheelpos(2, find(data.trialIDW == maxFramesIDW));
      reshapedData.timeW = 0:1/data.fps:max(timeWraw);
      maxFrames = length(reshapedData.timeW);
      reshapedData.wheelP = nan(maxFrames, Ntrials);
      reshapedData.wheelMovements = nan(maxFrames, Ntrials);
      reshapedData.wheelMovementsAcc = nan(maxFrames, Ntrials);
      reshapedData.wheelMovementsSz = nan(maxFrames, Ntrials);
      
      
      % Assignment
      for it2 = 1:length(trialID)
        curTrial = trialID(it2);
        try
          % Pupil
          validFrames = find(data.trialIDP == curTrial);
          if(length(validFrames) <= 2)
            continue;
          end
          validTimes = find(~isnan(data.parea(validFrames)));
          reshapedData.pupilArea(:, it2) = interp1(data.timeP(validFrames(validTimes)), data.parea(validFrames(validTimes)), reshapedData.timeP, 'pchip', NaN);
          validTimes = find(~isnan(data.ppos(1, validFrames)));
          reshapedData.pupilX(:, it2) = interp1(data.timeP(validFrames(validTimes)), data.ppos(1, validFrames(validTimes)), reshapedData.timeP, 'pchip', NaN);  
          validTimes = find(~isnan(data.ppos(2, validFrames)));
          reshapedData.pupilY(:, it2) = interp1(data.timeP(validFrames(validTimes)), data.ppos(2, validFrames(validTimes)), reshapedData.timeP, 'pchip', NaN);
          % Saccades
          validFrames = find(data.trialIDS == curTrial);
          for it1 = 1:length(validFrames)
            [~, cFrame] = min(abs(reshapedData.timeP-data.pmvton(validFrames(it1))));
            reshapedData.saccades(cFrame, it2) = 1;
          end

          % DFF
          validFrames = find(data.trialIDF == curTrial);
          if(length(validFrames) <= 2)
            continue;
          end
          if(~interpMovie)
            reshapedData.DFF(1:length(validFrames), :, it2) = data.resp(:, validFrames)';
          else
            validTimes = find(~isnan(sum(data.resp(:, validFrames), 1)));
            reshapedData.DFF(:, :, it2) = interp1(data.timeF(validFrames(validTimes)), data.resp(:, validFrames(validTimes))', reshapedData.timeF, 'pchip', NaN);
          end
          % Now the lowpass
          if(~isempty(lpFilt))
            validFrames = find(~isnan(sum(reshapedData.DFF(:, :, it2),2)));
            if(length(validFrames) <= 2)
              continue;
            end
            reshapedData.DFF(validFrames, :, it2) = filtfilt(lpFilt, reshapedData.DFF(validFrames, :, it2));
          end

          % Wheel
          validFrames = find(data.trialIDW == curTrial);
          if(length(validFrames) <= 2)
            continue;
          end
          validTimes = find(~isnan(data.wheelpos(1, validFrames)));
          reshapedData.wheelP(:, it2) = interp1(data.wheelpos(2, validFrames(validTimes)), data.wheelpos(1, validFrames(validTimes)), reshapedData.timeW, 'pchip', NaN);
          validFrames = find(data.trialIDWM == curTrial);
          for it1 = 1:length(validFrames)
            [~, cFrame] = min(abs(reshapedData.timeW-data.wheelmvton(validFrames(it1))));
            reshapedData.wheelMovements(cFrame, it2) = 1;
            reshapedData.wheelMovementsAcc(cFrame, it2) = data.wheelmvtacc(validFrames(it1));
            reshapedData.wheelMovementsSz(cFrame, it2) = data.wheelmvtsz(validFrames(it1));
          end
        catch ME
          slog.warning('Something went wrong on trial: %d. Setting its data to NaN', curTrial);
          warning(ME.message)
          fnames = {'pupilArea', 'pupilX', 'pupilY', 'saccades', 'wheelP', 'wheelMovements', 'wheelMovementsAcc', 'wheelMovementsSz'};
          for it1 = 1:length(fnames)
            reshapedData.(fnames{it1})(:, it2) = NaN;
          end
          reshapedData.DFF(:, :, it2) = NaN;
        end
      end
      % Now equilibrate num of frames
      if(~isempty(params.maxTime))
        maxFrames = min([length(reshapedData.timeP), length(reshapedData.timeF), length(reshapedData.timeW), find(reshapedData.timeF < params.maxTime, 1, 'last')]);
      else
        maxFrames = min([length(reshapedData.timeP), length(reshapedData.timeF), length(reshapedData.timeW)]);
      end
      fnames = {'pupilArea', 'pupilX', 'pupilY', 'saccades', 'wheelP', 'wheelMovements', 'wheelMovementsAcc', 'wheelMovementsSz'};
      for it1 = 1:length(fnames)
        reshapedData.(fnames{it1}) = reshapedData.(fnames{it1})(1:maxFrames, :);
      end
      reshapedData.DFF = reshapedData.DFF(1:maxFrames, :, :);
      reshapedData.time = reshapedData.timeF(1:maxFrames);
      reshapedData = rmfield(reshapedData, 'timeW');
      reshapedData = rmfield(reshapedData, 'timeP');
      reshapedData = rmfield(reshapedData, 'timeF');
      % Let's save to file also
      save(fullfile(obj.rootFolder, obj.name, obj.dataFolder, sprintf('%s_sessionRaw_%d.mat', obj.name, sessionNumber)), '-struct' ,'data');
      save(fullfile(obj.rootFolder, obj.name, obj.dataFolder, sprintf('%s_session_%d.mat', obj.name, sessionNumber)), '-struct', 'reshapedData');
      data = reshapedData;
      if(params.bar)
        ncbar.close();
      end
    end
    %%
    function cleanedData = cleanupSessionTrials(obj, sessionNumber, varargin)
    % CLEANUPSESSIONTRIALS Cleans trials from a session - removing data with nans - missalignments - lengths - simulated saccades
    % USAGE:
    %   data = obj.pullSessionTrials(sessionNumber);
    %
    % INPUT arguments:
    %   sessionNumber - index of the session
    %
    % INPUT optional arguments ('key' followed by its value):
    %   removeFnans   - to remove initial nans in the DFF signals
    %   realignPrestim - to realign prestim times - entry value will be the start of the OL
    %   lastTimePoint  - to fix the time vector
    %   removeSimSaccades - to remove simulated saccades
    %  maxTime
    %
    % OUTPUT arguments:
    %   cleanedData - data from the trials after cleaning
    %
    % EXAMPLE:
    %   dset = wf.dataset.load(pwd, 'animal1');
    %   data = dset.pullSessionTrials(21);
      
      params.removeFnans = true;
      params.realignPrestim = []; %
      params.lastTimePoint = [];
      params.removeSimSaccades = true;
      params.simSaccadesTimeThreshold = 3;
      params.openLoopDuration = 1.5;
      params.verbose = true;
      
      params = parse_pv_pairs(params, varargin);
      
      if(exist(fullfile(obj.rootFolder, obj.name, obj.dataFolder, sprintf('%s_session_%d.mat', obj.name, sessionNumber)), 'file'))
        data = load(fullfile(obj.rootFolder, obj.name, obj.dataFolder, sprintf('%s_session_%d.mat', obj.name, sessionNumber)));
        slog('Session %d data succesfully loaded from file: %s', sessionNumber, fullfile(obj.rootFolder, obj.name, obj.dataFolder, sprintf('%s_session_%d.mat', obj.name, sessionNumber)));
      else
        slog('Session %d data not found on: %s', sessionNumber, fullfile(obj.rootFolder, obj.name, obj.dataFolder, sprintf('%s_session_%d.mat', obj.name, sessionNumber)));
        return;
      end
      
      Ntrials = size(data.DFF, 3);
      invalidTrial = false(Ntrials, 1);
      fNames = {'pupilArea', 'pupilX', 'pupilY', 'saccades', 'wheelMovements', 'wheelMovementsAcc', 'wheelMovementsSz', 'wheelP'};
   
      if(params.removeSimSaccades)
        simSaccades = find(obj.data.eventTimes.OpenLoopEndedTime(data.trialID)' >= params.simSaccadesTimeThreshold);
        for it0 = 1:length(simSaccades)
          slog('Trial: %d contains sim saccade', data.trialID(simSaccades(it0)));
          [~, lastGoodFrame] = min(abs(data.time'-(params.openLoopDuration+obj.data.eventTimes.OpenLoopStartedTime(data.trialID(simSaccades(it0))))));
          if(lastGoodFrame < length(data.time))
            for it1 = 1:length(fNames)
              data.(fNames{it1})((lastGoodFrame+1):end, simSaccades(it0)) = NaN;
            end
            data.DFF((lastGoodFrame+1):end, :, simSaccades(it0)) = NaN;
          end
        end
      end
      
      if(~isempty(params.realignPrestim))
        targetTime = params.realignPrestim;
        [~, targetFrame] = min(abs(data.time-targetTime));
        data.openLoopStartFrame = targetFrame;
        [~, tmpFrame] = min(abs(data.time'-obj.data.eventTimes.OpenLoopStartedTime(data.trialID)'));
        frameDist = targetFrame-tmpFrame;
        tmpData = data;

        % Permute the DFF
        for it0 = 1:length(frameDist)
          tmpData.DFF(:, :, it0) = circshift(tmpData.DFF(:, :, it0), frameDist(it0), 1);
          if(frameDist(it0) < 0)
            tmpData.DFF(end+(frameDist(it0)+1):end, :, it0) = NaN;
          elseif(frameDist(it0) > 0)
            tmpData.DFF(1:frameDist(it0), :, it0) = NaN;
          end
        end
        % Permute everything else
        for it1 = 1:length(fNames)
          for it0 = 1:length(frameDist)
            tmpData.(fNames{it1})(:, it0) = circshift(tmpData.(fNames{it1})(:, it0), frameDist(it0), 1);
            if(frameDist(it0) < 0)
              tmpData.(fNames{it1})(end+(frameDist(it0)+1):end, it0) = NaN;
            elseif(frameDist(it0) > 0)
              tmpData.(fNames{it1})(1:frameDist(it0), it0) = NaN;
            end
          end
        end
        data = tmpData;
      end
      
      if(params.removeFnans)
        for it1 = 1:Ntrials
          avgTrial = mean(data.DFF(:, :, it1), 2);
          firstNaN = find(isnan(avgTrial), 1, 'first');
          if(isempty(firstNaN))
            continue;
          else
            firstNotNaN = find(~isnan(avgTrial), 1, 'first');
            if(isempty(firstNotNaN) || firstNotNaN > firstNaN)
              invalidTrial(it1) = true;
              slog('Trial: %d marked for deletion (early NaNs)', data.trialID(it1));
            end
          end
        end
        % = {'pupilArea', 'pupilX', 'pupilY', 'saccades', 'wheelMovements', 'wheelMovementsAcc', 'wheelMovementsSz', 'wheelP'};
        for it0 = 1:length(fNames)
          data.(fNames{it0})(:, invalidTrial) = [];
        end
        data.DFF(:, :, invalidTrial) = [];
        data.trialID(invalidTrial) = [];
      end
      
      if(~isempty(params.lastTimePoint))
        dt = mean(diff(data.time));
        if(data.time(end)+dt <= params.lastTimePoint)
          newFrames = length(data.time(1):dt:params.lastTimePoint)-length(data.time);
          slog('Need to add %d frames to the session block', newFrames);
          data.time = data.time(1):dt:params.lastTimePoint;
          for it0 = 1:length(fNames)
            data.(fNames{it0}) = [data.(fNames{it0}); nan(newFrames, length(data.trialID))];
          end
          data.DFF = [data.DFF; nan(newFrames, size(data.DFF, 2), size(data.DFF, 3))];
        elseif(data.time(end) >= params.lastTimePoint)
          oldFrames = find(data.time > params.lastTimePoint);
          slog('Need to remove %d frames from the session block', length(oldFrames));
          data.time(oldFrames) = [];
          for it0 = 1:length(fNames)
            data.(fNames{it0})(oldFrames, :) = [];
          end
          data.DFF(oldFrames, :, :) = [];
        end
      end
      % Do a last pass on array sizes
      for it0 = 1:length(fNames)
        if(size(data.(fNames{it0}), 1) < length(data.time))
          slog('Need to add frames to %s', fNames{it0});
          data.(fNames{it0}) = [data.(fNames{it0}); nan(length(data.time)-size(data.(fNames{it0}), 1), length(data.trialID))];
        elseif(size(data.(fNames{it0}), 1) > length(data.time))
          slog('Need to remove frames from %s', fNames{it0});
          data.(fNames{it0})((length(data.time)+1):end,:) = [];
        end
      end
      % Now for DFF
      if(size(data.DFF, 1) < length(data.time))
        slog('Need to add frames to DFF');
        data.DFF = [data.DFF; nan(length(data.time)-size(data.DFF, 1), size(data.DFF, 2), length(data.trialID))];
      elseif(size(data.DFF, 1) > length(data.time))
        slog('Need to remove frames from DFF');
        data.DFF((length(data.time)+1):end, :, :) = [];
      end
      
      cleanedData = data;
      save(fullfile(obj.rootFolder, obj.name, obj.dataFolder, sprintf('%s_sessionClean_%d.mat', obj.name, sessionNumber)), '-struct' ,'cleanedData');
    end
    %%
    function fullData = joinCleanSessions(obj, varargin)
    % JOINCLEANSESSSIONS Joins data from clean sessions together
    % USAGE:
    %   data = obj.joinCleanSessions();
    %
    % INPUT optional arguments ('key' followed by its value):
    %   sessList - index of the sessions to link - if empty - will use all
    %   valid sessions
    %   verbose   - true/false
    %   bar       - true/false
    %
    % OUTPUT arguments:
    %   data - data from the joint
    %
    % EXAMPLE:
    %   dset = wf.dataset.load(pwd, 'animal1');
    %   data = dset.joinCleanSessions();
      
      params.sessList = [];
      params.verbose = true;
      params.bar = true;
      params = parse_pv_pairs(params, varargin);
      
      if(isempty(params.sessList))
        params.sessList = obj.getValidSessions();
        if(params.verbose)
          slog('sessList empty. Concatenating all valid sessions');
        end
      end
      sessList = params.sessList;
      if(params.bar)
        ncbar('Joining sessions');
      end
      for it2 = 1:length(sessList)
        data = obj.pullSessionTrials(sessList(it2), 'type', 'clean');
        % To remember where the trial came from on the concatenated set (even
        % tho this info is on the trialSummary too)
        data.sessionID = ones(size(data.trialID))*sessList(it2);
        if(it2 == 1)
          fullData = data;
          fNames = fieldnames(fullData);
        else
          % Concatenate the trials
          for it3 = 1:length(fNames)
            % Don't concatenate the time
            if(strcmp(fNames{it3}, 'time'))
              continue;
            end
            fullData.(fNames{it3}) = cat(ndims(fullData.(fNames{it3})), fullData.(fNames{it3}), data.(fNames{it3}));
          end
        end
        if(params.bar)
          ncbar.update(it2/length(sessList));
        end
      end
      % Make sure time is a row vector
      if(size(fullData.time, 2) > 1)
        fullData.time = fullData.time';
      end
      save(fullfile(obj.rootFolder, obj.name, obj.dataFolder, sprintf('%s_jointSessions.mat', obj.name)), '-struct' ,'fullData');
      if(params.bar)
        ncbar.close();
      end
    end
    %%
    function fullData = loadJointSessions(obj)
      fullData = load(fullfile(obj.rootFolder, obj.name, obj.dataFolder, sprintf('%s_jointSessions.mat', obj.name)));
    end
    %%
    function fullData = loadCompressedSessions(obj)
      fullData = load(fullfile(obj.rootFolder, obj.name, obj.dataFolder, sprintf('%s_jointSessionsCompressed.mat', obj.name)));
    end
    %%
    function baseImg = getLocaNMFcomponent(obj, locaNMF, compIndex)
      refFrame = imref2d(size(obj.data.areasAllenLocaNMF));
      mask = obj.data.areasPixelMask;
      maskRegistered  = imwarp(mask, obj.data.alignment.transform, 'nearest', 'OutputView', refFrame);
      maskValid = find(maskRegistered);
      
      baseImg = nan(size(maskRegistered));
      validPixels = maskValid;
      baseImg(validPixels) = 0;
      validPixels = (baseImg' == 0);

      baseImg = nan(size(maskRegistered))';
      if(compIndex == 0)
        validPixels = maskValid;
        baseImg(validPixels) = 0;
        validPixels = (baseImg == 0);
        baseImg(validPixels) = locaNMF.means;
        baseImg = baseImg';
      else
        baseImg(validPixels) = locaNMF.compXY(compIndex, :)';
      end
      baseImg = baseImg';
    end
    %%
    function registeredData = registerCompressedSessions(obj, referenceImg, varargin)
    % REGISTERCOMPRESSEDSESSIONS Registers the data to the new ref
    % USAGE:
    %   data = obj.joinCleanSessions();
    %
    % INPUT optional arguments ('key' followed by its value):
    %   bar       - true/false
    %
    % OUTPUT arguments:
    %   data - data from the joint
    %
    % EXAMPLE:
    %   dset = wf.dataset.load(pwd, 'animal1');
    %   data = dset.registerCompressedSessions();
      params.bar = true;
      params = parse_pv_pairs(params, varargin);
      
      if(params.bar)
         ncbar.automatic('Loading compressed data');
      end
      fullData = obj.loadCompressedSessions();
      if(params.bar)
        ncbar.unsetAutomaticBar();
        ncbar.setBarTitle('Registering frames');
      end

      registeredRef = zeros(referenceImg.ImageSize);
      registeredData = struct;
      registeredData.DFF = nan(size(fullData.DFF, 1), numel(registeredRef));
      registeredData.frames = fullData.frames;
      registeredData.trialID = fullData.trialID;
      registeredData.time = fullData.time;
      registeredData.reference = referenceImg;
      registeredData.transform = obj.data.alignment.transform;
      originalSize = size(obj.data.areas);
      
      for it1 = 1:size(fullData.DFF, 1)
        newImg = imwarp(reshape(fullData.DFF(it1, :), originalSize), registeredData.transform, 'OutputView', registeredData.reference);
        registeredData.DFF(it1, :) = newImg(:);
        if(params.bar)
          ncbar.update(it1/size(fullData.DFF, 1));
        end
      end
      if(params.bar)
        ncbar.setBarTitle('Saving registereddata');
        ncbar.setAutomaticBar();
      end
      save(fullfile(obj.rootFolder, obj.name, obj.dataFolder, sprintf('%s_jointSessionsRegistered.mat', obj.name)), '-struct' ,'registeredData');
      if(params.bar)
        ncbar.close();
      end
    end
    %%
    function fullData = loadRegisteredSessions(obj)
      fullData = load(fullfile(obj.rootFolder, obj.name, obj.dataFolder, sprintf('%s_jointSessionsRegistered.mat', obj.name)));
    end
    %%
    function data = SVDcontinuousSessions(obj, varargin)
      params.bar = true;
      params.verbose = true;
      params.targetVar = [];
      params.export = 'pdf';
      params.mainTag = 'pcaVarExplained';
      params.appendTag = '';
      params = parse_pv_pairs(params, varargin);
      
      % Generate the mask of valid pixels
      refFrame = imref2d(size(obj.data.areasAllenLocaNMF));
      mask = obj.data.areasPixelMask;
      maskRegistered  = imwarp(mask, obj.data.alignment.transform, 'nearest', 'OutputView', refFrame);
      maskValid = find(maskRegistered);
      % Load the data
      if(params.bar)
         ncbar.automatic('Loading registered data');
      end
      data = obj.loadRegisteredSessions();
      data.DFF = data.DFF(:, maskValid); % Only care about these pixels
      if(params.bar)
        ncbar.setBarTitle('Running SVD');
      end
      
      %  Transpose and remove to save memory
      DFFt = data.DFF'; % Now it's pixels x Time
      data = rmfield(data, 'DFF');
      
      % Now PCA through SVD
      tabulatedValues = 0.05:0.05:1; % These come from https://arxiv.org/abs/1305.5870
      tabulatedCts = [1.5066 1.5816  1.6466  1.7048 1.7580 1.8074 1.8537 1.8974 1.9389 1.9786 2.0167 2.0533 2.0887 2.1229 2.1561 2.1883 2.2197 2.2503 2.2802 2.3094];

      % 0-mean data
      meanList = mean(DFFt, 2);
      DFFt = DFFt-meanList;
     
      if(params.verbose)
        slog('Running SVD');
      end
      
      [V, S, U] = svd(DFFt, 'econ'); % Appears to use less memory
      if(params.verbose)
        slog('Done!');
      end
      
      eigenv = diag(S).^2;
      latent = eigenv/sum(eigenv);
      % The components plot
      cumLatent = cumsum(latent);
      try
        ct = interp1(tabulatedValues, tabulatedCts, size(DFFt,1)/size(DFFt,2), 'pchip');
        largestComponent = find(latent <= mean(latent)*ct, 1, 'first');
        hFig = figure;
        plot(cumLatent/sum(latent),'.-');
        hold on;
        xl = xlim;
        h1 = plot(xl, [cumLatent(largestComponent); 0.95; 0.99]*[1 1], '--', 'LineWidth', 2);

        legend(h1, {sprintf('%.0f%%: %d comps (noise floor)', cumLatent(largestComponent)*100, largestComponent), ...
          sprintf('95%%: %d comps', find(cumsum(latent)/sum(latent) >= 0.95, 1, 'first')), ...
          sprintf('99%%: %d comps', find(cumsum(latent)/sum(latent) >= 0.99, 1, 'first'))}, 'location', 'SE', 'box', 'off');
        xlabel('# components');
        ylabel('Variance explained');
        offsetAxes(gca, 50);
        title('Variance explained per # comps');
        
        obj.exportFig(hFig, params);
      catch ME
        slog.warning('Something went wrong while plotting variance explained');
        slog.warning(getReport(ME, 'extended', 'hyperlinks', 'off'));
      end
      
      if(isempty(params.targetVar))
        ct = interp1(tabulatedValues, tabulatedCts, size(DFFt,1)/size(DFFt,2), 'pchip');
        largestComponent = find(latent <= mean(latent)*ct, 1, 'first');
        targetVar = cumLatent(largestComponent)*100;
      else
        targetVar = params.targetVar;
        largestComponent = find(cumLatent >= targetVar/100, 1, 'first');
      end
      
      coeff = V(:, 1:largestComponent);
      score = U(:, 1:largestComponent)*S(1:largestComponent, 1:largestComponent);
      fprintf('Keeping the first %d components (%.0f%% variance)\n', largestComponent, targetVar);

      data.coeff = coeff; % Old Uc
      data.score = score; % Old Vc
      % Such that Vc'*Uc' = real data (transposed),  score*coeff'+means recovers real data
      data.means = meanList;
      data.varExplained = params.targetVar;
      data.largestComponent = largestComponent;
      data.latent = latent;
      data.mask = maskRegistered;
      data.maskValid = maskValid;
      save(fullfile(obj.rootFolder, obj.name, obj.dataFolder, sprintf('%s_SVD.mat', obj.name)), '-struct' ,'data');
      if(params.bar)
        ncbar.close();
      end
    end
    %%
    function data = SVDregisteredSessions(obj, varargin)
      params.bar = true;
      params.verbose = true;
      params.targetVar = [];
      params.export = 'pdf';
      params.mainTag = 'pcaVarExplained';
      params.appendTag = '';
      params = parse_pv_pairs(params, varargin);
      
      % Generate the mask of valid pixels
      refFrame = imref2d(size(obj.data.areasAllenLocaNMF));
      mask = obj.data.areasPixelMask;
      maskRegistered  = imwarp(mask, obj.data.alignment.transform, 'nearest', 'OutputView', refFrame);
      maskValid = find(maskRegistered);
      % Load the data
      if(params.bar)
         ncbar.automatic('Loading registered data');
      end
      data = obj.loadRegisteredSessions();
      data.DFF = data.DFF(:, maskValid); % Only care about these pixels
      if(params.bar)
        ncbar.setBarTitle('Running SVD');
      end
      
      %  Transpose and remove to save memory
      DFFt = data.DFF'; % Now it's pixels x Time
      data = rmfield(data, 'DFF');
      
      % Now PCA through SVD
      tabulatedValues = 0.05:0.05:1; % These come from https://arxiv.org/abs/1305.5870
      tabulatedCts = [1.5066 1.5816  1.6466  1.7048 1.7580 1.8074 1.8537 1.8974 1.9389 1.9786 2.0167 2.0533 2.0887 2.1229 2.1561 2.1883 2.2197 2.2503 2.2802 2.3094];

      % 0-mean data
      meanList = mean(DFFt, 2);
      DFFt = DFFt-meanList;
     
      if(params.verbose)
        slog('Running SVD');
      end
      
      [V, S, U] = svd(DFFt, 'econ'); % Appears to use less memory
      if(params.verbose)
        slog('Done!');
      end
      
      eigenv = diag(S).^2;
      latent = eigenv/sum(eigenv);
      % The components plot
      cumLatent = cumsum(latent);
      try
        ct = interp1(tabulatedValues, tabulatedCts, size(DFFt,1)/size(DFFt,2), 'pchip');
        largestComponent = find(latent <= mean(latent)*ct, 1, 'first');
        hFig = figure;
        plot(cumLatent/sum(latent),'.-');
        hold on;
        xl = xlim;
        h1 = plot(xl, [cumLatent(largestComponent); 0.95; 0.99]*[1 1], '--', 'LineWidth', 2);

        legend(h1, {sprintf('%.0f%%: %d comps (noise floor)', cumLatent(largestComponent)*100, largestComponent), ...
          sprintf('95%%: %d comps', find(cumsum(latent)/sum(latent) >= 0.95, 1, 'first')), ...
          sprintf('99%%: %d comps', find(cumsum(latent)/sum(latent) >= 0.99, 1, 'first'))}, 'location', 'SE', 'box', 'off');
        xlabel('# components');
        ylabel('Variance explained');
        offsetAxes(gca, 50);
        title('Variance explained per # comps');
        
        obj.exportFig(hFig, params);
      catch ME
        slog.warning('Something went wrong while plotting variance explained');
        slog.warning(getReport(ME, 'extended', 'hyperlinks', 'off'));
      end
      
      if(isempty(params.targetVar))
        ct = interp1(tabulatedValues, tabulatedCts, size(DFFt,1)/size(DFFt,2), 'pchip');
        largestComponent = find(latent <= mean(latent)*ct, 1, 'first');
        targetVar = cumLatent(largestComponent)*100;
      else
        targetVar = params.targetVar;
        largestComponent = find(cumLatent >= targetVar/100, 1, 'first');
      end
      
      coeff = V(:, 1:largestComponent);
      score = U(:, 1:largestComponent)*S(1:largestComponent, 1:largestComponent);
      fprintf('Keeping the first %d components (%.0f%% variance)\n', largestComponent, targetVar);

      data.coeff = coeff; % Old Uc
      data.score = score; % Old Vc
      % Such that Vc'*Uc' = real data (transposed),  score*coeff'+means recovers real data
      data.means = meanList;
      data.varExplained = params.targetVar;
      data.largestComponent = largestComponent;
      data.latent = latent;
      data.mask = maskRegistered;
      data.maskValid = maskValid;
      save(fullfile(obj.rootFolder, obj.name, obj.dataFolder, sprintf('%s_SVD.mat', obj.name)), '-struct' ,'data');
      if(params.bar)
        ncbar.close();
      end
    end
    %%
    function fullData = loadBehaviorSessions(obj)
      fullData = load(fullfile(obj.rootFolder, obj.name, obj.dataFolder, sprintf('%s_jointSessionsBehavior.mat', obj.name)));
    end
    %%
    function fullData = splitAndCompressJointSessions(obj, varargin)
      params.bar = true;
      params = parse_pv_pairs(params, varargin);
      if(params.bar)
        ncbar.automatic('Loading full data');
      end
      fullData = load(fullfile(obj.rootFolder, obj.name, obj.dataFolder, sprintf('%s_jointSessions.mat', obj.name)));
      if(params.bar)
        ncbar.unsetAutomaticBar();
        ncbar.setBarTitle('Compressing DF/F');
      end
      avgF = squeeze(nanmean(fullData.DFF, 2));
      totalFrames = sum(~isnan(avgF(:)));
      compressedData = struct;
      compressedData.DFF = nan(totalFrames, size(fullData.DFF, 2));
      compressedData.frames = nan(totalFrames, 1);
      compressedData.trialID = nan(totalFrames, 1);
      compressedData.time = fullData.time;
      curFrame = 1;
      for it1 = 1:size(avgF, 2)
        cFrames = find(~isnan(avgF(:, it1)));
        compressedData.DFF(curFrame:(curFrame+length(cFrames)-1), :) = fullData.DFF(cFrames, :, it1);
        compressedData.frames(curFrame:(curFrame+length(cFrames)-1)) = cFrames;
        compressedData.trialID(curFrame:(curFrame+length(cFrames)-1)) = fullData.trialID(it1);
        curFrame = curFrame+length(cFrames);
        if(params.bar)
          ncbar.update(it1/size(avgF, 2));
        end
      end
      fullData = rmfield(fullData, 'DFF');
      if(params.bar)
        ncbar.setBarTitle('Saving compressed data');
        ncbar.setAutomaticBar();
      end
      save(fullfile(obj.rootFolder, obj.name, obj.dataFolder, sprintf('%s_jointSessionsCompressed.mat', obj.name)), '-struct' ,'compressedData');
      save(fullfile(obj.rootFolder, obj.name, obj.dataFolder, sprintf('%s_jointSessionsBehavior.mat', obj.name)), '-struct' ,'fullData');
      if(params.bar)
        ncbar.close();
      end
    end
    %%
    function [fullData, t] = loadTimeSeries(obj, timeRange, varargin)
    % LOADTIMESERIES Loads existing time series for a given time Range
    % USAGE:
    %   data = obj.loadTimeSeries(trialNumber);
    %
    % INPUT arguments:
    %   timeRange     - time [start finish] of time series to load
    %
    % INPUT optional arguments ('key' followed by its value):
    %   mask          - mask to use to average the time series (default empty. 'areas')
    %   handles       - file handles to use (sessions) on the load. Default: empty (meaning all)
    %   useAnimalTime - If true, will use the animal time instead as reference (decision time). Then timeRange is considered a single value of time before the decision to use. (Default: false)
    %
    % OUTPUT arguments:
    %   fullData      - the time series
    %   t             - the time vector
    %
    % EXAMPLE:
    %   dset = wf.dataset.load(pwd, 'animal1');
    %   dataset = dataset.generateTimeSeries(0, 1, 2.5, 12.5);
    %   data = dataset.loadTimeSeries([0 2]);
      
      params.mask = [];
      params.handles = [];
      params.useAnimalTime = false;
      params = parse_pv_pairs(params, varargin);
      
      if(isempty(timeRange))
        timeRange = [obj.timeSeriesInfo.timePrestimulus(1) obj.timeSeriesInfo.timeClosedLoop(end)];
      end
      if(params.useAnimalTime)
        usingPrestimulus = true;
        usingOpenLoop = true;
        usingClosedLoop = true;
      else
        usingPrestimulus = false;
        usingOpenLoop = false;
        usingClosedLoop = false;
        if(any(obj.timeSeriesInfo.timePrestimulus >= timeRange(1) & obj.timeSeriesInfo.timePrestimulus < timeRange(2)))
          usingPrestimulus = true;
        end
        if(any(obj.timeSeriesInfo.timeOpenLoop >= timeRange(1) & obj.timeSeriesInfo.timeOpenLoop < timeRange(2)))
          usingOpenLoop = true;
        end
        if(any(obj.timeSeriesInfo.timeClosedLoop >= timeRange(1) & obj.timeSeriesInfo.timeClosedLoop < timeRange(2)))
          usingClosedLoop = true;
        end
      end
      dataPreStimulus = struct;
      dataPreStimulus.pX = [];
      dataPreStimulus.pY = [];
      dataPreStimulus.pA = [];
      dataPreStimulus.pS = [];
      dataPreStimulus.wP = [];
      dataPreStimulus.wV = [];
      dataPreStimulus.T = [];
      dataPreStimulus.neural = [];
      dataOpenLoop = dataPreStimulus;
      dataClosedLoop = dataPreStimulus;
      data = dataPreStimulus;
      data.index = [];
      data.decisionTime = [];
      emptyData = data;
      fullData = data;
      fnames = fieldnames(dataPreStimulus);
      if(~isempty(params.mask))
        switch params.mask
          case 'areas'
            % Mask.map contains the image with the coordinates of each area
            % Mask.list is a cell list with the pixel list of each area
            mask = struct;
            mask.map = obj.data.trialsSummary.aux.visarea_ed;
            mapCount = hist(mask.map(:), 1:max(mask.map(:)));
            invalid = find(mapCount < 2);
            for it = 1:length(invalid)
              fprintf('Warning, removing area %d because of its size: %d.\n', invalid(it), mapCount(invalid(it)));
              mask.map(mask.map == invalid(it)) = 0;
            end
            mask.list = arrayfun(@(x)find(mask.map == x), unique(mask.map(~~mask.map)), 'UniformOutput', false);
            %tmpel = cellfun(@(x)numel(x), mask.list);
          otherwise
            mask = [];
        end
      else
        mask = [];
      end
      if(isempty(params.handles))
        validHandles = 1:size(obj.timeSeriesInfo.handles, 1);
      else
        validHandles = params.handles;
      end
      timeVector = [];
      if(usingPrestimulus)
        timeVector = [timeVector, obj.timeSeriesInfo.timePrestimulus];
      end
      if(usingOpenLoop)
        timeVector = [timeVector, obj.timeSeriesInfo.timeOpenLoop];
      end
      if(usingClosedLoop)
        timeVector = [timeVector, obj.timeSeriesInfo.timeClosedLoop];
      end
      if(params.useAnimalTime)
        timeVectorOriginal = timeVector;
        timeVector = -timeVector(end:-1:1);
        validFrames = find(timeVector >= -timeRange);
        t = timeVector(validFrames);
      else
        validFrames = find(timeVector >= timeRange(1) & timeVector < timeRange(2));
        t = timeVector(validFrames);
      end
      data.decisionTime = obj.timeSeriesInfo.decisionTime(data.index);
      data.decisionTimeCorr = obj.timeSeriesInfo.decisionTime(data.index);
      
      ncbar('Loading time series');
      for it = validHandles
        data = emptyData;
        data.index = find(obj.timeSeriesInfo.handleFile == it);
        if(usingPrestimulus)
          tmp = load(fullfile(obj.outputFolder, obj.timeSeriesInfo.handles{it, 1}));
          if(isstruct(mask))
            neural = cellfun(@(x)mean(tmp.neural(x, :, :), 1), mask.list, 'UniformOutput', false);
            tmp.neural = cell2mat(neural);
          end
          for it2 = 1:length(fnames)
            dataPreStimulus.(fnames{it2}) = tmp.(fnames{it2});
          end
        end
        if(usingOpenLoop)
          tmp = load(fullfile(obj.outputFolder, obj.timeSeriesInfo.handles{it, 2}));
          if(isstruct(mask))
            neural = cellfun(@(x)mean(tmp.neural(x, :, :)), mask.list, 'UniformOutput', false);
            tmp.neural = cell2mat(neural);
          end
          for it2 = 1:length(fnames)
            dataOpenLoop.(fnames{it2}) = tmp.(fnames{it2});
          end
        end
        if(usingClosedLoop)
          tmp = load(fullfile(obj.outputFolder, obj.timeSeriesInfo.handles{it, 3}));
          if(isstruct(mask))
            neural = cellfun(@(x)mean(tmp.neural(x, :, :)), mask.list, 'UniformOutput', false);
            tmp.neural = cell2mat(neural);
          end
          for it2 = 1:length(fnames)
            try
            dataClosedLoop.(fnames{it2}) = tmp.(fnames{it2});
            catch
              it
              fnames{it2}
              size(dataClosedLoop.(fnames{it2}))
              size(tmp.(fnames{it2}))
            end
          end
        end
        
        % Concatenate all data
        for it2 = 1:length(fnames)
          if(usingPrestimulus)
            switch fnames{it2}
              case 'neural'
                data.(fnames{it2}) = cat(2, data.(fnames{it2}), dataPreStimulus.(fnames{it2}));
              otherwise
                data.(fnames{it2}) = cat(1, data.(fnames{it2}), dataPreStimulus.(fnames{it2}));
            end
          end
          if(usingOpenLoop)
            switch fnames{it2}
              case 'neural'
                data.(fnames{it2}) = cat(2, data.(fnames{it2}), dataOpenLoop.(fnames{it2}));
              otherwise
                data.(fnames{it2}) = cat(1, data.(fnames{it2}), dataOpenLoop.(fnames{it2}));
            end
          end
          if(usingClosedLoop)
            switch fnames{it2}
              case 'neural'
                data.(fnames{it2}) = cat(2, data.(fnames{it2}), dataClosedLoop.(fnames{it2}));
              otherwise
                data.(fnames{it2}) = cat(1, data.(fnames{it2}), dataClosedLoop.(fnames{it2}));
            end
          end
        end
        % Fix decision time for animal time
        if(params.useAnimalTime)
          for it2 = 1:length(data.index)
            try
              posTime = find(~isnan(sum(data.neural(:, :, it2),1)) & sum(data.neural(:, :, it2),1) ~= 0, 1, 'last');
              if(~isempty(posTime))
                data.decisionTimeCorr(it2) = timeVectorOriginal(posTime);
              else
                data.decisionTimeCorr(it2) = NaN;
              end
            catch ME
              data.decisionTimeCorr(it2) = NaN;
              fprintf('%s\n', ME.message);
            end
          end
        end
        % Trim time series
        for it2 = 1:length(fnames)
          if(params.useAnimalTime)
            switch fnames{it2}
              case 'neural'
                newData = nan(size(data.(fnames{it2}),1), length(validFrames), size(data.(fnames{it2}),3));
                for it3 = 1:length(data.index)
                  curValidFrames = find(timeVectorOriginal >= data.decisionTimeCorr(it3)-timeRange & timeVectorOriginal < data.decisionTimeCorr(it3));
                  newData(:, (end-length(curValidFrames)+1):end, it3) = data.(fnames{it2})(:, curValidFrames, it3);
                end
                data.(fnames{it2}) = newData;
              otherwise
                newData = nan(length(validFrames), size(data.(fnames{it2}),2));
                for it3 = 1:length(data.index)
                  curValidFrames = find(timeVectorOriginal >= data.decisionTimeCorr(it3)-timeRange & timeVectorOriginal < data.decisionTimeCorr(it3));
                  newData((end-length(curValidFrames)+1):end, it3) = data.(fnames{it2})(curValidFrames, it3);
                end
                data.(fnames{it2}) = newData;
            end
          else
            switch fnames{it2}
              case 'neural'
                data.(fnames{it2}) = data.(fnames{it2})(:, validFrames, :);
              otherwise
                data.(fnames{it2}) = data.(fnames{it2})(validFrames, :);
            end
          end
        end
        % Now concatenate across sessions
        for it2 = 1:length(fnames)
          switch fnames{it2}
            case 'neural'
              fullData.(fnames{it2}) = cat(ndims(data.(fnames{it2})), fullData.(fnames{it2}), data.(fnames{it2}));
            otherwise
              fullData.(fnames{it2}) = cat(ndims(data.(fnames{it2})), fullData.(fnames{it2}), data.(fnames{it2}));
          end
        end
        fullData.index = [fullData.index; data.index];
        
        ncbar.update(find(it==validHandles)/length(validHandles));
        % Things to concatenate
        %data.index
        %dataPreStimulus.(fnames{it2})
 %       dataPreStimulus.(fnames{it2}) = cat(ndims(tmp.(fnames{it2})), ...
%                                                dataPreStimulus.(fnames{it2}), tmp.(fnames{it2}));
%   dataOpenLoop.(fnames{it2}) = cat(ndims(tmp.(fnames{it2})), ...
%                                              dataOpenLoop.(fnames{it2}), tmp.(fnames{it2}));
%             dataClosedLoop.(fnames{it2}) = cat(ndims(tmp.(fnames{it2})), ...
%                                                dataClosedLoop.(fnames{it2}), tmp.(fnames{it2}));
      end

      ncbar.close();
    end
    
  end
end
