% CD to this folder
addpath(genpath(fullfile(pwd,'helperFunctions')));

%% PANEL A
clearvars;
load('data/fig5a.mat');

hFig = createCenteredFigure('width', 21, 'height', 5);
p = panel();
p.margin = 12;
p.pack(1, [0.24 -1]);

p(1, 1).select();

cmap = cblindmap();
hold on;
errorbar(angleList, ylowm, ylowl, ylowu, '-o', 'MarkerEdgeColor', 'k', 'MarkerFaceColor', cmap(1, :), 'Color', cmap(1, :), 'MarkerSize', 4)
errorbar(angleList, yhighm, yhighl, yhighu', '-o', 'MarkerEdgeColor', 'k', 'MarkerFaceColor', cmap(2, :), 'Color', cmap(2, :), 'MarkerSize', 4)
l = legend('low attention', 'high attention');
l.Location = 'southeast';
l.Box = 'off';
l.Position(1) = l.Position(1) + 0.05;
set(gca, 'XTick', -90:30:90);
set(gca, 'YTick', 0:0.2:1);
spaceOutAxes();
offsetAxes(gca);
xlabel('Angle difference |\theta_L|-|\theta_R|');
ylabel('Prob of choosing right');
set(findall(hFig,'-property','FontName'),'FontName', 'Arial');
set(findall(hFig,'-property','FontSize'),'FontSize', 8);

%% PANELS B C D
clearvars;
load('data/fig5bcd.mat');
projLabels = {'left', 'right'};

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
plot([1 1], yl, 'k--');
xlabel('time (s)');
ylabel('Projection (a.u.)');
xlim([0.5 2.5]);
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
plot([1 1], yl, 'k--');
xlabel('time (s)');
ylabel('d''');
xlim([0.5 2.5]);
offsetAxes();
spaceOutAxes();

% d' spatial map
p(1, 4).select();

cmap = parula(256);
x = imagesc(CC, 'XData', [0 2.5], 'YData', [0 2.5]);
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
set(gca, 'XTick', 0:2);
set(gca, 'YTick', 0:2);
plot([1 1], yl, 'k--');
plot(xl, [1 1], 'k--');
xlabel('Time from mov (s)');
ylabel('Time from mov (s)');

set(findall(hFig,'-property','FontName'),'FontName', 'Arial');
set(findall(hFig,'-property','FontSize'),'FontSize', 8);

%% PANELS E F

clearvars;
load('data/fig5e.mat');
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

load('data/fig5f.mat');
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

%% PANEL G
clearvars;
load('data/fig5g.mat');

vectorModes = {'left E', 'right E', 'left H', 'right H'};
hFig = createCenteredFigure('width', 21, 'height', 7);
p = panel();
p.margin = 10;
p.pack(1, [0.05 0.5 0.35 0.2]);

%%% The projection
p(1 ,2).select();

hold on;

h = [];
for it0 = 1:8
  if(it0 < 5)
    h = [h; plot3(xyz{it0}(:,1), xyz{it0}(:,2), xyz{it0}(:,3), '-', 'Color', cidx{it0}, 'DisplayName', labels{it0})];
  else
    h = [h; plot3(xyz{it0}(:,1), xyz{it0}(:,2), xyz{it0}(:,3), '--', 'Color', cidx{it0}, 'DisplayName', labels{it0})];
  end
end
  
grid on;
set(gca, 'XTick', -0.5:0.5:2);
xlim([-0.5 1.5]);
set(gca, 'YTick', -1:1:1);
view([26.4 34.1]);
xlabel('time from stimulus (s)');
ylabel('choice');
zlabel('attention');
l = legend(h, 'Location', 'eastoutside');
l.Box = 'off';
set(findall(hFig,'-property','FontName'),'FontName', 'Arial');
set(findall(hFig,'-property','FontSize'),'FontSize', 8);

%% PANEL H
clearvars;
load('data/fig5h.mat');
hFig = createCenteredFigure('width', 21, 'height', 5);
p = panel();
p.margin = 12;
p.pack(1, [0.24 -1]);
p(1,1).select();

imagesc(diffList, diffList, diffMat);
set(gca, 'XTick', 0:30:90);
set(gca, 'YTick', 0:30:90);
xlabel('Difficulty ||\theta_L-\theta_R||');
ylabel('Difficulty ||\theta_L-\theta_R||');
axis equal square tight ij;
cb = colorbar;
caxis([0 30]);
cb.Label.String = 'Angle between choice vectors(°)';
xlim([0 90]);
ylim([0 90]);

set(findall(hFig,'-property','FontName'),'FontName', 'Arial');
set(findall(hFig,'-property','FontSize'),'FontSize', 8);

%% PANEL I

clearvars;
load('data/fig5i.mat');
hFig = createCenteredFigure('width', 21, 'height', 5);
p = panel();
p.margin = 12;
p.pack(1, [0.24 -1]);
p(1,1).select();
imagesc(attList, attList, diffMat);
set(gca, 'XTick', 0:0.25:1);
set(gca, 'YTick', 0:0.25:1);
xlabel('Attention level');
ylabel('Attention level');
axis equal square tight ij;
cb = colorbar;
caxis([0 30]);
cb.Label.String = 'Angle between choice vectors(°)';
xlim([0 1]);
ylim([0 1]);

set(findall(hFig,'-property','FontName'),'FontName', 'Arial');
set(findall(hFig,'-property','FontSize'),'FontSize', 8);

%% PANEL J
clearvars;
load('data/fig5j.mat');

cmap = parula(length(attList)+1);
hFig = createCenteredFigure('width', 21, 'height', 5);
p = panel();
p.margin = 12;
p.pack(1, [0.24 -1]);
p(1,1).select();
hold on;
for it1 = 1:length(attList)
  plot(diffList, attListPsy(:, it1), '-', 'MarkerEdgeColor', 'k', 'MarkerFaceColor', cmap(it1, :), 'Color', cmap(it1, :))
end

set(gca, 'XTick', -90:30:90);
set(gca, 'YTick', 0:0.2:1);
spaceOutAxes();
offsetAxes(gca);
xlabel('Angle difference |\theta_L|-|\theta_R|');
ylabel('Prob of choosing right');
cb = colorbar;
cb.Label.String = 'Attention';
set(findall(hFig,'-property','FontName'),'FontName', 'Arial');
set(findall(hFig,'-property','FontSize'),'FontSize', 8);