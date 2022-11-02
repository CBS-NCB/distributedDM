%%% -----------------------------------------------------------------------
%%% This script will generate panels from Figure 1
%%% -----------------------------------------------------------------------
%%% Requires `main.m` to be run first.

clearvars;
addpath(genpath(pwd));
projectFolder = [pwd filesep '..' filesep 'widefieldChoice' filesep 'data' filesep 'excMice'];
load([projectFolder filesep 'projectData.mat']);

%% Psychometric curves Right choice

angleDivisions = 13;
fullCurve = [];
hFig = createCenteredFigure('width',5,'height', 4);
hold on;
cmap = lines(1);
for itAnimal = 1:Ndatasets
  dset = wf.dataset.load(projectFolder, datasetList.dataset{itAnimal});
  dset.rootFolder = [pwd filesep 'data' filesep 'excMice']; % Overwrite for machine swaps
  
  sessList = dset.getValidSessions();
  validTrials = ismember(dset.data.trialsSummary.sess, sessList);
  choiceDir = dset.data.trialsSummary.choicedir(validTrials);
  outcome = dset.data.trialsSummary.outcome(validTrials);

  validAngles = -diff(abs(dset.data.trialsSummary.ori), [],2); % Left - Right
  angleList = unique(validAngles);
  [~, binnedAngles, angleListIdx] = histcounts(angleList, linspace(-90, 90, angleDivisions));
  binnedAnglesCenter = binnedAngles(1:end-1)+diff(binnedAngles)/2;

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
    %validAngles = -diff(abs(dset.data.trialsSummary.ori(validTrials, :)), [],2); % Left - Right
    if(datasetList.reference{itAnimal} == 90)
      validAngles = diff(abs(dset.data.trialsSummary.ori(validTrials, :)), [],2); % Left - Right
    else
      validAngles = -diff(abs(dset.data.trialsSummary.ori(validTrials, :)), [],2); % Left - Right
    end
    %validAngles = -diff(abs(dset.data.trialsSummary.ori(validTrials, :)), [],2); % Left - Right
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

  [~, pci] = binofit(round(rightChoice.*rightChoiceN), rightChoiceN);
  sem = diff(pci,[],2)/2;
  m = pci(:,1)+sem;
  h = plot(binnedAnglesCenter, m*100, '-', 'Color', [cmap(1, :) 0.5]);
  fullCurve = [fullCurve, m];
  h.LineWidth = 0.5;
end
  
m = mean(fullCurve, 2);
sem = std(fullCurve, [], 2)/sqrt(Ndatasets);
h = errorbar(binnedAnglesCenter, m*100,sem*100, '-o', 'MarkerSize', 4, 'DisplayName', 'global', 'Color', cmap(1, :));
h.MarkerFaceColor =  'k';


xlabel('Angle difference \langle|\theta_L|-|\theta_R|\rangle');
ylabel('Right Choice (%)');
xlim([-90 90]);
ylim([0 100]);
set(gca, 'XTick', linspace(-90, 90, 5));
set(gca, 'YTick', linspace(0, 100,5));
spaceOutAxes(gca);
offsetAxes(gca, 50);
set(findall(hFig,'-property','FontName'),'FontName', 'Arial');
set(findall(hFig,'-property','FontSize'),'FontSize', 8);
set(hFig,'Color','w');
exportgraphics(hFig, 'psychRight.pdf', 'ContentType', 'vector');

%% Psychometric curves time outs

angleDivisions = 13;
fullCurve = [];
hFig = createCenteredFigure('width',5,'height', 4);
hold on;
cmap = lines(1);
for itAnimal = 1:Ndatasets
  dset = wf.dataset.load(projectFolder, datasetList.dataset{itAnimal});
  dset.rootFolder = [pwd filesep 'data' filesep 'excMice']; % Overwrite for machine swaps
  % Assign the session list
 %dset.plotPsychometricCurves();
  sessList = dset.getValidSessions();
  validTrials = ismember(dset.data.trialsSummary.sess, sessList);
  choiceDir = dset.data.trialsSummary.choicedir(validTrials);
  outcome = dset.data.trialsSummary.outcome(validTrials);

  validAngles = -diff(abs(dset.data.trialsSummary.ori), [],2); % Left - Right
  angleList = unique(validAngles);
  [~, binnedAngles, angleListIdx] = histcounts(angleList, linspace(-90, 90, angleDivisions));
  binnedAnglesCenter = binnedAngles(1:end-1)+diff(binnedAngles)/2;

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
    if(datasetList.reference{itAnimal} == 90)
      validAngles = diff(abs(dset.data.trialsSummary.ori(validTrials, :)), [],2); % Left - Right
    else
      validAngles = -diff(abs(dset.data.trialsSummary.ori(validTrials, :)), [],2); % Left - Right
    end
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

  [~, pci] = binofit(round(timeOut.*timeOutN), timeOutN);
  sem = diff(pci,[],2)/2;
  m = pci(:,1)+sem;
  h = plot(binnedAnglesCenter, m*100, '-', 'Color', [cmap(1, :) 0.5]);
  fullCurve = [fullCurve, m];
  h.LineWidth = 0.5;
end

  
m = mean(fullCurve, 2);
sem = std(fullCurve, [], 2)/sqrt(Ndatasets);
h = errorbar(binnedAnglesCenter, m*100,sem*100, '-o', 'MarkerSize', 4, 'DisplayName', 'global', 'Color', cmap(1, :));
h.MarkerFaceColor =  'k';


xlabel('Angle difference \langle|\theta_L|-|\theta_R|\rangle');
ylabel('Time out (%)');
xlim([-90 90]);
ylim([0 20]);
set(gca, 'XTick', linspace(-90, 90, 5));
set(gca, 'YTick', linspace(0, 20,5));
spaceOutAxes(gca);
offsetAxes(gca, 50);
set(findall(hFig,'-property','FontName'),'FontName', 'Arial');
set(findall(hFig,'-property','FontSize'),'FontSize', 8);
set(hFig,'Color','w');
exportgraphics(hFig, 'psychTimeout.pdf', 'ContentType', 'vector');

