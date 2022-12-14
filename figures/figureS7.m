% CD to this folder
addpath(genpath(fullfile(pwd,'helperFunctions')));

%%
clearvars;
load('data/figs7.mat');

hFig = createCenteredFigure('width', 21, 'height', 6);
p = panel();
p.margin = 10;
p.pack(1, [0.28 0.28 0.16 0.28]);

p(1, 4).select();

hold on;
% Plot parameters
compressionLevel = 60;
BarWidth = 1;
WhiskersWidthRatio = 8;
i = 1;
j = 2;
subCmap = lines(2);
XPositions = 1:2;
lineStyle = 'all';
pointSize = 16;
pointOpacity = 0.7;
pointBorderSize = 0.1;
fullData = d_data;
for k = [i, j]
yValues = fullData(:, k);
yMean = nanmean(yValues);
yStd = nanstd(yValues);
ySem = yStd/sqrt(size(yValues,1));
yCI = ySem*1.96;
%plot the standard deviation box
 p1 = patch('XData', [0 1 1 0]*2*BarWidth/WhiskersWidthRatio+XPositions(k)-BarWidth/WhiskersWidthRatio, ...
    'YData', [0 0 1 1]*2*yStd+yMean-yStd, 'FaceColor',subCmap(i,:),'EdgeColor', subCmap(i,:),'LineWidth',0.1);
%plot the mean+-SEM box
p2 = patch('XData', [0 1 1 0]*2*BarWidth/WhiskersWidthRatio+XPositions(k)-BarWidth/WhiskersWidthRatio, ...
      'YData', [0 0 1 1]*2*yCI+yMean-yCI, 'FaceColor',subCmap(i,:),'EdgeColor', subCmap(i,:),'LineWidth',0.1);
%plot the mean line 
% Manual alpha-like feature
alphaColor1 = subCmap(i,:)*0.2+0.8;
alphaColor1(alphaColor1 > 1) = 1;
alphaColor2 = subCmap(i,:)*0.35+0.65;
alphaColor2(alphaColor2 > 1) = 1;
set(p1, 'FaceColor', alphaColor1, 'EdgeColor', alphaColor1);
set(p2, 'FaceColor', alphaColor2, 'EdgeColor', alphaColor2);
plot([XPositions(k)-BarWidth/WhiskersWidthRatio XPositions(k)+BarWidth/WhiskersWidthRatio],[yMean yMean],'Color',subCmap(i,:),'LineWidth',2)
end

plot([XPositions(i) XPositions(j)],[nanmean(fullData(:,i)),nanmean(fullData(:,j))],'Color',subCmap(i,:), 'LineWidth',1)

% The actual line and scatter
plot(repmat([XPositions(i), XPositions(j)], size(fullData, 1), 1)', [fullData(:, i), fullData(:, j)]', 'Color', subCmap(i, :));

hl = scatter(XPositions(i)*ones(size(fullData,1), 1), fullData(:, i), pointSize, 'o','MarkerEdgeColor', subCmap(i,:)*pointOpacity, 'MarkerFaceColor',subCmap(i,:),'LineWidth',pointBorderSize);
scatter(XPositions(j)*ones(size(fullData,1), 1), fullData(:, j), pointSize, 'o','MarkerEdgeColor', subCmap(i,:)*pointOpacity, 'MarkerFaceColor',subCmap(i,:),'LineWidth',pointBorderSize);
[a,b]=ttest(fullData(1,:)-fullData(2,:));
hh = sigstar([1 2], b);
hh{1}(2).Position(2) = hh{1}(2).Position(2)*1.01;
xlim([0 3]);
set(gca,'XTick', 1:2);
set(gca,'XTickLabel', {sprintf('OL movement'), 'CL movement'});
ylabel('Average performance');
ylim([0.6 0.8]);
offsetAxes();
spaceOutAxes();


p(1, 3).select();
hold on;
% Plot parameters
compressionLevel = 60;
BarWidth = 1;
WhiskersWidthRatio = 8;
i = 1;
subCmap = lines(2);
XPositions = 1:2;
lineStyle = 'all';
pointSize = 16;
pointOpacity = 0.7;
pointBorderSize = 0.1;
fullData = c_data;
 for k = i
  yValues = fullData(:, k);
  yMean = nanmean(yValues);
  yStd = nanstd(yValues);
  ySem = yStd/sqrt(size(yValues,1));
  yCI = ySem*1.96;
  %plot the standard deviation box
  p1 = patch('XData', [0 1 1 0]*2*BarWidth/WhiskersWidthRatio+XPositions(k)-BarWidth/WhiskersWidthRatio, ...
      'YData', [0 0 1 1]*2*yStd+yMean-yStd, 'FaceColor',subCmap(i,:),'EdgeColor', subCmap(i,:),'LineWidth',0.1);
      %'YData', [0 0 1 1]*2*yStd+yMean-yStd, 'FaceColor',StdColor(StdIndex,:),'EdgeColor', StdColor(StdIndex,:),'LineWidth',0.1);
  %plot the mean+-SEM box
  p2 = patch('XData', [0 1 1 0]*2*BarWidth/WhiskersWidthRatio+XPositions(k)-BarWidth/WhiskersWidthRatio, ...
        'YData', [0 0 1 1]*2*yCI+yMean-yCI, 'FaceColor',subCmap(i,:),'EdgeColor', subCmap(i,:),'LineWidth',0.1);
  %plot the mean line 
  % Manual alpha-like feature
  alphaColor1 = subCmap(i,:)*0.2+0.8;
  alphaColor1(alphaColor1 > 1) = 1;
  alphaColor2 = subCmap(i,:)*0.35+0.65;
  alphaColor2(alphaColor2 > 1) = 1;
  set(p1, 'FaceColor', alphaColor1, 'EdgeColor', alphaColor1);
  set(p2, 'FaceColor', alphaColor2, 'EdgeColor', alphaColor2);
  plot([XPositions(k)-BarWidth/WhiskersWidthRatio XPositions(k)+BarWidth/WhiskersWidthRatio],[yMean yMean],'Color',subCmap(i,:),'LineWidth',2)
end    

hl = scatter(XPositions(i)*ones(size(fullData,1), 1), fullData(:, i), pointSize, 'o','MarkerEdgeColor', subCmap(i,:)*pointOpacity, 'MarkerFaceColor',subCmap(i,:),'LineWidth',pointBorderSize);

xlim([0 2]);
set(gca,'XTick', 1);
set(gca,'XTickLabel', 'Same movement');
%set(gca,'XTickLabel', {sprintf('OL movement'), 'CL movement'});
ylabel('Percentage (%%)');
%ylim([0.6 0.8]);
offsetAxes();
spaceOutAxes();

p(1,1).select();
hold on;


i = 1;
j = 2;
fullData = a_data;
for k = [i, j]
  yValues = fullData(:, k);
  yMean = nanmean(yValues);
  yStd = nanstd(yValues);
  ySem = yStd/sqrt(size(yValues,1));
  yCI = ySem*1.96;
  %plot the standard deviation box
   p1 = patch('XData', [0 1 1 0]*2*BarWidth/WhiskersWidthRatio+XPositions(k)-BarWidth/WhiskersWidthRatio, ...
      'YData', [0 0 1 1]*2*yStd+yMean-yStd, 'FaceColor',subCmap(i,:),'EdgeColor', subCmap(i,:),'LineWidth',0.1);
  %plot the mean+-SEM box
  p2 = patch('XData', [0 1 1 0]*2*BarWidth/WhiskersWidthRatio+XPositions(k)-BarWidth/WhiskersWidthRatio, ...
        'YData', [0 0 1 1]*2*yCI+yMean-yCI, 'FaceColor',subCmap(i,:),'EdgeColor', subCmap(i,:),'LineWidth',0.1);
  %plot the mean line 
  % Manual alpha-like feature
  alphaColor1 = subCmap(i,:)*0.2+0.8;
  alphaColor1(alphaColor1 > 1) = 1;
  alphaColor2 = subCmap(i,:)*0.35+0.65;
  alphaColor2(alphaColor2 > 1) = 1;
  set(p1, 'FaceColor', alphaColor1, 'EdgeColor', alphaColor1);
  set(p2, 'FaceColor', alphaColor2, 'EdgeColor', alphaColor2);
  plot([XPositions(k)-BarWidth/WhiskersWidthRatio XPositions(k)+BarWidth/WhiskersWidthRatio],[yMean yMean],'Color',subCmap(i,:),'LineWidth',2)
end    

plot([XPositions(i) XPositions(j)],[nanmean(fullData(:,i)),nanmean(fullData(:,j))],'Color',subCmap(i,:), 'LineWidth',1)
plot(repmat([XPositions(i), XPositions(j)], size(fullData, 1), 1)', [fullData(:, i), fullData(:, j)]', 'Color', subCmap(i, :));

hl = scatter(XPositions(i)*ones(size(fullData,1), 1), fullData(:, i), pointSize, 'o','MarkerEdgeColor', subCmap(i,:)*pointOpacity, 'MarkerFaceColor',subCmap(i,:),'LineWidth',pointBorderSize);
scatter(XPositions(j)*ones(size(fullData,1), 1), fullData(:, j), pointSize, 'o','MarkerEdgeColor', subCmap(i,:)*pointOpacity, 'MarkerFaceColor',subCmap(i,:),'LineWidth',pointBorderSize);
[a,b]=ttest(fullData(:,1)-fullData(:,2));
hh = sigstar([1 2], b);

hh{1}(2).Position(2) = hh{1}(2).Position(2)*1.01;
xlim([0 3]);
set(gca,'XTick', 1:2);
set(gca,'XTickLabel', {'low attention', 'high attention'});
ylabel('Reaction Time (s)');
offsetAxes();
spaceOutAxes();


p(1,2).select();

hold on;

fullData = b_data;
i = 1;
j = 2;
for k = [i, j]
  yValues = fullData(:, k);
  yMean = nanmean(yValues);
  yStd = nanstd(yValues);
  ySem = yStd/sqrt(size(yValues,1));
  yCI = ySem*1.96;
  %plot the standard deviation box
   p1 = patch('XData', [0 1 1 0]*2*BarWidth/WhiskersWidthRatio+XPositions(k)-BarWidth/WhiskersWidthRatio, ...
      'YData', [0 0 1 1]*2*yStd+yMean-yStd, 'FaceColor',subCmap(i,:),'EdgeColor', subCmap(i,:),'LineWidth',0.1);
  %plot the mean+-SEM box
  p2 = patch('XData', [0 1 1 0]*2*BarWidth/WhiskersWidthRatio+XPositions(k)-BarWidth/WhiskersWidthRatio, ...
        'YData', [0 0 1 1]*2*yCI+yMean-yCI, 'FaceColor',subCmap(i,:),'EdgeColor', subCmap(i,:),'LineWidth',0.1);
  %plot the mean line 
  % Manual alpha-like feature
  alphaColor1 = subCmap(i,:)*0.2+0.8;
  alphaColor1(alphaColor1 > 1) = 1;
  alphaColor2 = subCmap(i,:)*0.35+0.65;
  alphaColor2(alphaColor2 > 1) = 1;
  set(p1, 'FaceColor', alphaColor1, 'EdgeColor', alphaColor1);
  set(p2, 'FaceColor', alphaColor2, 'EdgeColor', alphaColor2);
  plot([XPositions(k)-BarWidth/WhiskersWidthRatio XPositions(k)+BarWidth/WhiskersWidthRatio],[yMean yMean],'Color',subCmap(i,:),'LineWidth',2)
end    
plot([XPositions(i) XPositions(j)],[nanmean(fullData(:,i)),nanmean(fullData(:,j))],'Color',subCmap(i,:), 'LineWidth',1)
plot(repmat([XPositions(i), XPositions(j)], size(fullData, 1), 1)', [fullData(:, i), fullData(:, j)]', 'Color', subCmap(i, :));
hl = scatter(XPositions(i)*ones(size(fullData,1), 1), fullData(:, i), pointSize, 'o','MarkerEdgeColor', subCmap(i,:)*pointOpacity, 'MarkerFaceColor',subCmap(i,:),'LineWidth',pointBorderSize);
scatter(XPositions(j)*ones(size(fullData,1), 1), fullData(:, j), pointSize, 'o','MarkerEdgeColor', subCmap(i,:)*pointOpacity, 'MarkerFaceColor',subCmap(i,:),'LineWidth',pointBorderSize);
[a,b]=ttest(fullData(:,1)-fullData(:,2));
hh = sigstar([1 2], b);
hh{1}(2).Position(2) = hh{1}(2).Position(2)*1.01;
xlim([0 3]);
set(gca,'XTick', 1:2);
set(gca,'XTickLabel', {'low attention', 'high attention'});
ylabel('Average performance (%%)');
offsetAxes();
spaceOutAxes();

set(findall(hFig,'-property','FontName'),'FontName', 'Arial');
set(findall(hFig,'-property','FontSize'),'FontSize', 7);

