clearvars;
load('fig3hi.mat');

hFig = createCenteredFigure('width', 21, 'height', 27);
p = panel();
p.margin = 10;
p.pack([0.1 0.1 0.1 0.7], 1);
p(3,1).pack(1, [0.15 0.25 0.15 0.15]);

compressionLevel = 60;
BarWidth = 1;
WhiskersWidthRatio = 8;
XPositions = 1:2;

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
plot(repmat([XPositions(i), XPositions(j)], size(dataPairsH, 1), 1)', [dataPairsH(:, i), dataPairsH(:, j)]', 'Color', subCmap(3, :));
hl = scatter(XPositions(i)*ones(size(dataPairsH,1), 1), dataPairsH(:, i), pointSize, 'o','MarkerEdgeColor', subCmap(3,:)*pointOpacity, 'MarkerFaceColor',subCmap(3,:),'LineWidth',pointBorderSize);
scatter(XPositions(j)*ones(size(dataPairsH,1), 1), dataPairsH(:, j), pointSize, 'o','MarkerEdgeColor', subCmap(3,:)*pointOpacity, 'MarkerFaceColor',subCmap(3,:),'LineWidth',pointBorderSize);

[a,b]=ttest(dataPairsH(:,1)-dataPairsH(:,2));
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
set(gca, 'XTickLabel', {'low att.', 'high att.'});
xlim([0.5 2.5]);
ylim([0 2.1]);
ylabel('d''');
offsetAxes();
spaceOutAxes();



%%% The d prime
p(3,1,1,2).select();



hold on;
h1 = [];
h2 = [];
fill([t' fliplr(t')], [y1u' fliplr(y1l')], cmap(1,:), 'FaceAlpha', 0.5, 'EdgeColor', 'none');
h1 = [h1; plot(t, y1m, '-', 'Color', cmap(1, :), 'DisplayName', sprintf('%s', projLabels{1}))];
fill([t' fliplr(t')], [y2u' fliplr(y2l')], cmap(2,:), 'FaceAlpha', 0.5, 'EdgeColor', 'none');
h2 = [h2; plot(t, y2m, '-', 'Color', cmap(2, :), 'DisplayName', sprintf('%s', projLabels{2}))];

ylim([0 2]);
yl = ylim;
plot([0 0], yl, 'k--');
l = legend([h1 h2]);
l.Box = 'off';
l.Location = 'best';
xlabel('time (s)');

ylabel('d''');
xlim([-1 1]);
set(gca, 'XTick', -1:1);
offsetAxes();
spaceOutAxes();

% Now the d prime spatial map
p(3,1,1,3).select();

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
plot(repmat([XPositions(i), XPositions(j)], size(dataPairsI, 1), 1)', [dataPairsI(:, i), dataPairsI(:, j)]', 'Color', subCmap(3, :));
hl = scatter(XPositions(i)*ones(size(dataPairsI,1), 1), dataPairsI(:, i), pointSize, 'o','MarkerEdgeColor', subCmap(3,:)*pointOpacity, 'MarkerFaceColor',subCmap(3,:),'LineWidth',pointBorderSize);
scatter(XPositions(j)*ones(size(dataPairsI,1), 1), dataPairsI(:, j), pointSize, 'o','MarkerEdgeColor', subCmap(3,:)*pointOpacity, 'MarkerFaceColor',subCmap(3,:),'LineWidth',pointBorderSize);

[a,b]=ttest(dataPairsI(:,1)-dataPairsI(:,2));
if(b <= 0.05)
  hh = sigstar([1 2], b);
  hh{1}(2).Position(2) = hh{1}(2).Position(2)*1.01;
  hh{1}(1).YData = hh{1}(1).YData*1.01;
end
set(gca,'XTick', 1:2);
set(gca, 'XTickLabel', projLabels);
xlim([0.5 2.5]);
ylabel('d''');
title('High attention');
offsetAxes();
spaceOutAxes();

set(findall(hFig,'-property','FontName'),'FontName', 'Arial');
set(findall(hFig,'-property','FontSize'),'FontSize', 8);
