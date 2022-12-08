function hList = plotPairedData(fullData, varargin)

params.cmap = cblindmap();
params.BarWidth = 0.5;
params.WhiskersWidthRatio = 2;
params.pointBorderSize = 0.5;
params.pointSize = 8;
params.pointOpacity = 0.8;
params.lineStyle = 'all';
params = parse_pv_pairs(params, varargin);

subCmap = params.cmap;
BarWidth = params.BarWidth;
WhiskersWidthRatio = params.WhiskersWidthRatio;
pointBorderSize = params.pointBorderSize;
pointSize = params.pointSize;
pointOpacity = params.pointOpacity;
lineStyle = params.lineStyle;

XPositions = [1 2];
i = 1;
j = 2;
hold on;
for k = [i, j]
  yValues = fullData(:, k);
  yMean = nanmean(yValues);
  yStd = nanstd(yValues);
  ySem = yStd/sqrt(size(yValues,1));
  %plot the standard deviation box
   p1 = patch('XData', [0 1 1 0]*2*BarWidth/WhiskersWidthRatio+XPositions(k)-BarWidth/WhiskersWidthRatio, ...
      'YData', [0 0 1 1]*2*yStd+yMean-yStd, 'FaceColor',subCmap(i,:),'EdgeColor', subCmap(i,:),'LineWidth',0.1);
      %'YData', [0 0 1 1]*2*yStd+yMean-yStd, 'FaceColor',StdColor(StdIndex,:),'EdgeColor', StdColor(StdIndex,:),'LineWidth',0.1);
  %plot the mean line 
  % Manual alpha-like feature
  alphaColor1 = subCmap(i,:)*0.2+0.8;
  alphaColor1(alphaColor1 > 1) = 1;
  alphaColor2 = subCmap(i,:)*0.35+0.65;
  alphaColor2(alphaColor2 > 1) = 1;
  set(p1, 'FaceColor', alphaColor1, 'EdgeColor', alphaColor1);
  plot([XPositions(k)-BarWidth/WhiskersWidthRatio XPositions(k)+BarWidth/WhiskersWidthRatio],[yMean yMean],'Color',subCmap(i,:),'LineWidth',2)
end    
switch lineStyle
  case 'all'
    plot(repmat([XPositions(i), XPositions(j)], size(fullData, 1), 1)', [fullData(:, i), fullData(:, j)]', 'Color', subCmap(i, :));
  case 'mean'
    plot([XPositions(i) XPositions(j)],[nanmean(fullData(:,i)),nanmean(fullData(:,j))],'Color',subCmap(i,:), 'LineWidth',1)
end

hList = scatter(XPositions(i)*ones(size(fullData,1), 1), fullData(:, i), pointSize, 'o','MarkerEdgeColor', subCmap(i,:)*pointOpacity, 'MarkerFaceColor',subCmap(i,:),'LineWidth',pointBorderSize);
hList = scatter(XPositions(j)*ones(size(fullData,1), 1), fullData(:, j), pointSize, 'o','MarkerEdgeColor', subCmap(i,:)*pointOpacity, 'MarkerFaceColor',subCmap(i,:),'LineWidth',pointBorderSize);
