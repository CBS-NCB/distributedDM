% CD to this folder
addpath(genpath(fullfile(pwd,'helperFunctions')));

%% PANELS A B C
clearvars;
load('data/fig3abc.mat');

hFig = createCenteredFigure('width', 21, 'height', 5);
p = panel();
p.margin = 10;
p.pack(1, [0.1 0.35 0.35 0.2]);

%%% Projection trajs
p(1 ,2).select();
cmap = cblindmap();
hold on;
h1 = [];
h2 = [];
fill([t' fliplr(t')], [y1u' fliplr(y1l')], cmap(1,:), 'FaceAlpha', 0.5, 'EdgeColor', 'none');
h1 = [h1; plot(t, y1m, '-', 'Color', cmap(1, :), 'DisplayName', sprintf('%s', projLabels{1}))];
fill([t' fliplr(t')], [y2u' fliplr(y2l')], cmap(2,:), 'FaceAlpha', 0.5, 'EdgeColor', 'none');
h2 = [h2; plot(t, y2m, '-', 'Color', cmap(2, :), 'DisplayName', sprintf('%s', projLabels{2}))];
yl = ylim;
%ylim([0 yl(2)*1.2]);
plot([0 0], yl, 'k--');
xlabel('time (s)');
ylabel('Projection (a.u.)');
xlim([-1 1]);
l = legend([h1 h2]);
l.Box = 'off';
l.Location = 'best';
offsetAxes();
spaceOutAxes();

%%% d' traj
p(1, 3).select();
hold on;
fill([t' fliplr(t')], [ydu' fliplr(ydl')], cmap(3,:), 'FaceAlpha', 0.5, 'EdgeColor', 'none');
plot(t, ydm, '-', 'Color', cmap(3, :));
ylim([0 1.6])
set(gca, 'YTick', 0:0.5:1.5);
yl = ylim;
plot([0 0], yl, 'k--');
xlabel('time (s)');
ylabel('d''');
xlim([-1 1]);
offsetAxes();
spaceOutAxes();

% d' spatial map
p(1, 4).select();

cmap = parula(256);
x = imagesc(CC, 'XData', [-1 1], 'YData', [-1 1]);
x.AlphaData = ~isnan(CC);
setImageAxis();
colormap(cmap);
axp = get(gca, 'Position');
c = colorbar;
c.Label.String = 'Correlation';
caxis([0 1]);
axis square;
hold on;
xl = xlim;
yl = ylim;
plot([0 0], yl, 'k--');
plot(xl, [0 0], 'k--');
xlabel('Time from mov (s)');
ylabel('Time from mov (s)');

set(findall(hFig,'-property','FontName'),'FontName', 'Arial');
set(findall(hFig,'-property','FontSize'),'FontSize', 8);

%% PANELS D E F G
clearvars;
load('data/fig3defg.mat');

hFig = createCenteredFigure('width', 21, 'height', 27);
p = panel();
p.margin = 10;
p.pack([0.1 0.1 0.1 0.7], 1);

p(2,1).pack(1, [0.25 0.25 0.25 0.25]);
p(2,1,1,1).select();

cmapOrig = [0 0 0; cblindmap];
cmap = cmapOrig(regionColor, :);
hold on;
for it1 = 1:length(regionNames)
  h = plot(areaDt, areaD(:, it1), 'Color', cmap(it1, :));
end

xlabel('time from mov (s)');
ylabel('d''');
ylim([0 1.2]);
plot([0 0], ylim, 'k--');
spaceOutAxes(gca, 0.05);
offsetAxes(gca);

p(2,1,1,2).select();
hold on;
for it1 = 1:length(regionNames)
  h1 = plot(areaDfitX, areaDfit(:, it1), 'Color', cmap(it1, :), 'LineWidth', 0.5);
end

xlabel('time from mov (s)');
ylabel('d''');
xlim([-1 0.2]);
ylim([0 1.2]);
plot([0 0], ylim, 'k--');
spaceOutAxes(gca, 0.05);
offsetAxes(gca);

p(2,1,1,3).select();
hold on;
cmapOrig = [0 0 0; cblindmap];
cmap = cmapOrig([1 regionColor], :);
errorbar(slopeDx, slopeDm, slopeDl, slopeDu, 'o', 'Color', cmap(3,:), 'MarkerSize', 4);

xlim([0 max(regionPosition)+1]);
set(gca, 'XTick', [1 regionPosition]);
set(gca, 'XTickLabel', {'global', regionNames{regionOrder}});
set(gca, 'XTickLabelRotation', 90);
plot(xlim,[0 0], 'k--');
offsetAxes();
spaceOutAxes();
ylabel('Slope pre movement')

p(2,1,1,4).select();
hold on;
errorbar(slopeDx, slopeTm, slopeTl, slopeTu, 'o', 'Color', cmap(4,:), 'MarkerSize', 4);

xlim([0 max(regionPosition)+1]);
set(gca, 'XTick', [1 regionPosition]);
set(gca, 'XTickLabel', {'global', regionNames{regionOrder}});
set(gca, 'XTickLabelRotation', 90);
plot(xlim,[0 0], 'k--');
offsetAxes();

set(findall(hFig,'-property','FontName'),'FontName', 'Arial');
set(findall(hFig,'-property','FontSize'),'FontSize', 8);

%% PANELS H I
clearvars;
load('data/fig3hi.mat');

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

%% PANEL J
clearvars;
load('data/fig3j.mat');

hFig = createCenteredFigure('width', 5, 'height', 5);
cmap = cblindmap();
cmap = cmap(1:end,:);
hold on;
  
%%% The projection
for it0 = 1:4
  switch it0
    case 1
      attMult = 1;
      diffMult = 1;
    case 2
      attMult = 1;
      diffMult = 0;
    case 3
      attMult = 0;
      diffMult = 1;
    case 4
      attMult = 0;
      diffMult = 0;
  end
  opacity = (1+diffMult)/2;
  h1 = plot3(x(:, (it0-1)*2+1), y(:, (it0-1)*2+1), z(:, (it0-1)*2+1), '-', 'Color', cmap(1, :)*opacity);
  h2 = plot3(x(:, (it0-1)*2+2), y(:, (it0-1)*2+2), z(:, (it0-1)*2+2), '-', 'Color', cmap(1+1, :)*opacity);
end

yl = ylim;
plot([0 0], ylim, 'k--');
set(gca, 'XTick', 0:2);
xlabel('time (s)');
ylabel('Projection (a.u.)');
xlim([-0.5 0.5]);
set(gca,'XTick',-0.5:0.5:1);
set(gca,'YTick',-2:1:2);
ylim([-2 2]);
zlim([0 1.2]);
grid on;
set(gca,'ZTick',0:0.25:1);
view([34.5 30.5]);
