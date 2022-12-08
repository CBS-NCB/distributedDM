clearvars;
load('fig5e.mat');
labels = {'easy', 'hard'};
hFig = createCenteredFigure('width', 21, 'height', 27);
p = panel();
p.margin = 10;
p.pack([0.1 0.1 0.1 0.7], 1);
p(3,1).pack(1, [0.15 0.15 0.15 0.15]);

compressionLevel = 60;
BarWidth = 1;
WhiskersWidthRatio = 8;

lineStyle = 'all';
pointSize = 16;
pointOpacity = 0.7;
pointBorderSize = 0.1;
cmap = cblindmap();
cmap = cmap(1:end,:);
subCmap = cmap;

p(3,1,1,1).select();

i = 1;
j = 2;
XPositions = [1 2];
hold on;
for k = [i, j]

   p1 = patch('XData', [0 1 1 0]*2*BarWidth/WhiskersWidthRatio+XPositions(k)-BarWidth/WhiskersWidthRatio, ...
      'YData', [cil(k) cil(k) ciu(k) ciu(k)], 'FaceColor',subCmap(3,:),'EdgeColor', subCmap(3,:),'LineWidth',0.1);
  % Manual alpha-like feature
  alphaColor1 = subCmap(3,:)*0.2+0.8;
  alphaColor1(alphaColor1 > 1) = 1;
  alphaColor2 = subCmap(3,:)*0.35+0.65;
  alphaColor2(alphaColor2 > 1) = 1;
  set(p1, 'FaceColor', alphaColor1, 'EdgeColor', alphaColor1);
  plot([XPositions(k)-BarWidth/WhiskersWidthRatio XPositions(k)+BarWidth/WhiskersWidthRatio],[m(k) m(k)],'Color',subCmap(3,:),'LineWidth',2)
end    
plot([XPositions(i) XPositions(j)],[m(i) m(j)],'Color',subCmap(3,:), 'LineWidth',1)
plot(repmat([XPositions(i), XPositions(j)], size(dataPairs, 1), 1)', [dataPairs(:, i), dataPairs(:, j)]', 'Color', subCmap(3, :));
hl = scatter(XPositions(i)*ones(size(dataPairs,1), 1), dataPairs(:, i), pointSize, 'o','MarkerEdgeColor', subCmap(3,:)*pointOpacity, 'MarkerFaceColor',subCmap(3,:),'LineWidth',pointBorderSize);
scatter(XPositions(j)*ones(size(dataPairs,1), 1), dataPairs(:, j), pointSize, 'o','MarkerEdgeColor', subCmap(3,:)*pointOpacity, 'MarkerFaceColor',subCmap(3,:),'LineWidth',pointBorderSize);

[a,b]=ttest(dataPairs(:,1)-dataPairs(:,2));
if(b <= 0.05)
  hh = sigstar([1 2], b);
  hh{1}(2).Position(2) = hh{1}(2).Position(2)*1.01;
  hh{1}(1).YData = hh{1}(1).YData*1.01;
else
  hh = sigstar([1 2], b);
  hh{1}(2).Position(2) = hh{1}(2).Position(2)*1.15;
  hh{1}(1).YData = hh{1}(1).YData*1.1;
end
set(gca,'XTick', 1:2);
set(gca, 'XTickLabel', labels);
xlim([0.5 2.5]);
ylim([0 4]);
ylabel('d''');
title('High attention');
offsetAxes();
spaceOutAxes();

load('fig5f.mat');
labels = {'easy', 'hard'};

p(3,1,1,2).select();

i = 1;
j = 2;
XPositions = [1 2];
hold on;
for k = [i, j]

   p1 = patch('XData', [0 1 1 0]*2*BarWidth/WhiskersWidthRatio+XPositions(k)-BarWidth/WhiskersWidthRatio, ...
      'YData', [cil(k) cil(k) ciu(k) ciu(k)], 'FaceColor',subCmap(3,:),'EdgeColor', subCmap(3,:),'LineWidth',0.1);
  % Manual alpha-like feature
  alphaColor1 = subCmap(3,:)*0.2+0.8;
  alphaColor1(alphaColor1 > 1) = 1;
  alphaColor2 = subCmap(3,:)*0.35+0.65;
  alphaColor2(alphaColor2 > 1) = 1;
  set(p1, 'FaceColor', alphaColor1, 'EdgeColor', alphaColor1);
  plot([XPositions(k)-BarWidth/WhiskersWidthRatio XPositions(k)+BarWidth/WhiskersWidthRatio],[m(k) m(k)],'Color',subCmap(3,:),'LineWidth',2)
end    
plot([XPositions(i) XPositions(j)],[m(i) m(j)],'Color',subCmap(3,:), 'LineWidth',1)
plot(repmat([XPositions(i), XPositions(j)], size(dataPairs, 1), 1)', [dataPairs(:, i), dataPairs(:, j)]', 'Color', subCmap(3, :));
hl = scatter(XPositions(i)*ones(size(dataPairs,1), 1), dataPairs(:, i), pointSize, 'o','MarkerEdgeColor', subCmap(3,:)*pointOpacity, 'MarkerFaceColor',subCmap(3,:),'LineWidth',pointBorderSize);
scatter(XPositions(j)*ones(size(dataPairs,1), 1), dataPairs(:, j), pointSize, 'o','MarkerEdgeColor', subCmap(3,:)*pointOpacity, 'MarkerFaceColor',subCmap(3,:),'LineWidth',pointBorderSize);

[a,b]=ttest(dataPairs(:,1)-dataPairs(:,2));
if(b <= 0.05)
  hh = sigstar([1 2], b);
  hh{1}(2).Position(2) = hh{1}(2).Position(2)*1.01;
  hh{1}(1).YData = hh{1}(1).YData*1.01;
else
  hh = sigstar([1 2], b);
  hh{1}(2).Position(2) = hh{1}(2).Position(2)*1.15;
  hh{1}(1).YData = hh{1}(1).YData*1.1;
end
set(gca,'XTick', 1:2);
set(gca, 'XTickLabel', labels);
xlim([0.5 2.5]);
ylim([0 4]);
ylabel('d''');
title('Low attention');
offsetAxes();
spaceOutAxes();