% CD to this folder
addpath(genpath(fullfile(pwd,'helperFunctions')));

%% PANEL A
clearvars;
load('data/fig2a.mat');

%save('fig2a.mat', 'contourListTransform', 'labelCenterListTransform', 'areasNames', 'V1', 'RL', 'PM');

hFig = createCenteredFigure('width', 5, 'height', 5);
hold on;
for it3 = 1:length(contourListTransform)
  plot(contourListTransform{it3}(:,1),contourListTransform{it3}(:,2),'k');
end
setImageAxis();
xlim([1 67]);
ylim([1 67]);
box on;

h =[];
for it1 = 1:length(V1)
  ax = gca;
  ax.ColorOrderIndex = 1;
  for it2 = 1:3
    switch it2
      case 1
        plot(V1{it1}(:,1), V1{it1}(:,2));
      case 2
        plot(PM{it1}(:,1), PM{it1}(:,2));
      case 3
        plot(RL{it1}(:,1), RL{it1}(:,2));
    end
  end
end
for it3 = 1:length(contourListTransform)
  if(labelCenterListTransform{it3}(1) < 67 & labelCenterListTransform{it3}(2) < 67 & labelCenterListTransform{it3}(1) > 0 & labelCenterListTransform{it3}(2) > 0)
    text(labelCenterListTransform{it3}(1), labelCenterListTransform{it3}(2), areasNames{it3}, 'HorizontalAlignment', 'center');
  end
end

set(findall(hFig,'-property','FontName'),'FontName', 'Arial');
set(findall(hFig,'-property','FontSize'),'FontSize', 8);

%% PANEL D
clearvars;
load('data/fig2d.mat');

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
ylim([0 yl(2)*1.2]);
yl = ylim;
plot([0 0], yl, 'k--');
hp1 = fill([barX(1) barX(1) barX(end) barX(end)], [yl(2)*0.95 yl(2) yl(2) yl(2)*0.95], [1 1 1]*0.5, 'FaceAlpha', 0.5, 'EdgeColor', 'none');
hp1.FaceAlpha = 0.5;
uistack(hp1, 'bottom');
xlabel('time (s)');
ylabel('Projection (a.u.)');
xlim(xRange);
l = legend([h1 h2]);
l.Box = 'off';
l.Location = 'best';
offsetAxes();
spaceOutAxes();

%%% d' traj
p(1, 3).select();
hold on;
fill([td' fliplr(td')], [ydu' fliplr(ydl')], cmap(3,:), 'FaceAlpha', 0.5, 'EdgeColor', 'none');
plot(td, ydm, '-', 'Color', cmap(3, :));
ylim([0 1.6])
set(gca, 'YTick', 0:0.5:1.5);
yl = ylim;
plot([0 0], yl, 'k--');
hp1 = fill([barX(1) barX(1) barX(end) barX(end)], [yl(2)*0.95 yl(2) yl(2) yl(2)*0.95], [1 1 1]*0.5, 'FaceAlpha', 0.5, 'EdgeColor', 'none');
hp1.FaceAlpha = 0.5;
uistack(hp1, 'bottom');
xlabel('time (s)');
ylabel('d''');
xlim(xRange);
offsetAxes();
spaceOutAxes();

% d' spatial map
p(1, 4).select();

cmap = parula(256);
x = imagesc(movData);
x.AlphaData = ~isnan(movData);
setImageAxis();
colormap(cmap);
axp = get(gca, 'Position');
c = colorbar;
c.Label.String = 'd''';

caxis([0 1.2]);

hold on;
for it3 = 1:length(contourListTransform)
  plot(contourListTransform{it3}(:,1),contourListTransform{it3}(:,2),'w');
end
for it3 = 1:length(contourListTransform)
  if(labelCenterListTransform{it3}(1) < 67 & labelCenterListTransform{it3}(2) < 67 & labelCenterListTransform{it3}(1) > 0 & labelCenterListTransform{it3}(2) > 0)
    text(labelCenterListTransform{it3}(1), labelCenterListTransform{it3}(2), areasNames{it3}, 'HorizontalAlignment', 'center', 'Color', [0 0 0]);
  end
end
axis off;
xlim([1 67]);
ylim([1 67]);  
set(gca, 'Position', axp-[0.025 0 0 0]);

set(findall(hFig,'-property','FontName'),'FontName', 'Arial');
set(findall(hFig,'-property','FontSize'),'FontSize', 8);

%% PANEL E
clearvars;
load('data/fig2e.mat');

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
ylim([0 yl(2)*1.2]);
yl = ylim;
plot([0 0], yl, 'k--');
hp1 = fill([barX(1) barX(1) barX(end) barX(end)], [yl(2)*0.95 yl(2) yl(2) yl(2)*0.95], [1 1 1]*0.5, 'FaceAlpha', 0.5, 'EdgeColor', 'none');
hp1.FaceAlpha = 0.5;
uistack(hp1, 'bottom');
xlabel('time (s)');
ylabel('Projection (a.u.)');
xlim(xRange);
l = legend([h1 h2]);
l.Box = 'off';
l.Location = 'best';
offsetAxes();
spaceOutAxes();

%%% d' traj
p(1, 3).select();
hold on;
fill([td' fliplr(td')], [ydu' fliplr(ydl')], cmap(3,:), 'FaceAlpha', 0.5, 'EdgeColor', 'none');
plot(td, ydm, '-', 'Color', cmap(3, :));
ylim([0 1.6])
set(gca, 'YTick', 0:0.5:1.5);
yl = ylim;
plot([0 0], yl, 'k--');
hp1 = fill([barX(1) barX(1) barX(end) barX(end)], [yl(2)*0.95 yl(2) yl(2) yl(2)*0.95], [1 1 1]*0.5, 'FaceAlpha', 0.5, 'EdgeColor', 'none');
hp1.FaceAlpha = 0.5;
uistack(hp1, 'bottom');
xlabel('time (s)');
ylabel('d''');
xlim(xRange);
offsetAxes();
spaceOutAxes();

% d' spatial map
p(1, 4).select();

cmap = parula(256);
x = imagesc(movData);
x.AlphaData = ~isnan(movData);
setImageAxis();
colormap(cmap);
axp = get(gca, 'Position');
c = colorbar;
c.Label.String = 'd''';

caxis([0 1.2]);

hold on;
for it3 = 1:length(contourListTransform)
  plot(contourListTransform{it3}(:,1),contourListTransform{it3}(:,2),'w');
end
for it3 = 1:length(contourListTransform)
  if(labelCenterListTransform{it3}(1) < 67 & labelCenterListTransform{it3}(2) < 67 & labelCenterListTransform{it3}(1) > 0 & labelCenterListTransform{it3}(2) > 0)
    text(labelCenterListTransform{it3}(1), labelCenterListTransform{it3}(2), areasNames{it3}, 'HorizontalAlignment', 'center', 'Color', [0 0 0]);
  end
end
axis off;
xlim([1 67]);
ylim([1 67]);  
set(gca, 'Position', axp-[0.025 0 0 0]);

set(findall(hFig,'-property','FontName'),'FontName', 'Arial');
set(findall(hFig,'-property','FontSize'),'FontSize', 8);

%% PANEL F
clearvars;
load('data/fig2f.mat');


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
ylim([0 yl(2)*1.2]);
yl = ylim;
plot([0 0], yl, 'k--');
hp1 = fill([barX(1) barX(1) barX(end) barX(end)], [yl(2)*0.95 yl(2) yl(2) yl(2)*0.95], [1 1 1]*0.5, 'FaceAlpha', 0.5, 'EdgeColor', 'none');
hp1.FaceAlpha = 0.5;
uistack(hp1, 'bottom');
xlabel('time (s)');
ylabel('Projection (a.u.)');
xlim(xRange);
l = legend([h1 h2]);
l.Box = 'off';
l.Location = 'best';
offsetAxes();
spaceOutAxes();

%%% d' traj
p(1, 3).select();
hold on;
fill([td' fliplr(td')], [ydu' fliplr(ydl')], cmap(3,:), 'FaceAlpha', 0.5, 'EdgeColor', 'none');
plot(td, ydm, '-', 'Color', cmap(3, :));
ylim([0 1.6])
set(gca, 'YTick', 0:0.5:1.5);
yl = ylim;
plot([0 0], yl, 'k--');
hp1 = fill([barX(1) barX(1) barX(end) barX(end)], [yl(2)*0.95 yl(2) yl(2) yl(2)*0.95], [1 1 1]*0.5, 'FaceAlpha', 0.5, 'EdgeColor', 'none');
hp1.FaceAlpha = 0.5;
uistack(hp1, 'bottom');
xlabel('time (s)');
ylabel('d''');
xlim(xRange);
offsetAxes();
spaceOutAxes();

% d' spatial map
p(1, 4).select();

cmap = parula(256);
x = imagesc(movData);
x.AlphaData = ~isnan(movData);
setImageAxis();
colormap(cmap);
axp = get(gca, 'Position');
c = colorbar;
c.Label.String = 'd''';

caxis([0 1.2]);

hold on;
for it3 = 1:length(contourListTransform)
  plot(contourListTransform{it3}(:,1),contourListTransform{it3}(:,2),'w');
end
for it3 = 1:length(contourListTransform)
  if(labelCenterListTransform{it3}(1) < 67 & labelCenterListTransform{it3}(2) < 67 & labelCenterListTransform{it3}(1) > 0 & labelCenterListTransform{it3}(2) > 0)
    text(labelCenterListTransform{it3}(1), labelCenterListTransform{it3}(2), areasNames{it3}, 'HorizontalAlignment', 'center', 'Color', [0 0 0]);
  end
end
axis off;
xlim([1 67]);
ylim([1 67]);  
set(gca, 'Position', axp-[0.025 0 0 0]);

set(findall(hFig,'-property','FontName'),'FontName', 'Arial');
set(findall(hFig,'-property','FontSize'),'FontSize', 8);

%% PANEL G
clearvars;
load('data/fig2g.mat');

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
ylim([-1.2 yl(2)*1.2]);
yl = ylim;
plot([1 1], yl, 'k--');
hp1 = fill([barX(1) barX(1) barX(end) barX(end)], [yl(2)*0.95 yl(2) yl(2) yl(2)*0.95], [1 1 1]*0.5, 'FaceAlpha', 0.5, 'EdgeColor', 'none');
hp1.FaceAlpha = 0.5;
uistack(hp1, 'bottom');
xlabel('time (s)');
ylabel('Projection (a.u.)');
xlim(xRange);
l = legend([h1 h2]);
l.Box = 'off';
l.Location = 'best';
offsetAxes();
spaceOutAxes();

%%% d' traj
p(1, 3).select();
hold on;
fill([td' fliplr(td')], [ydu' fliplr(ydl')], cmap(3,:), 'FaceAlpha', 0.5, 'EdgeColor', 'none');
plot(td, ydm, '-', 'Color', cmap(3, :));
ylim([0 1.6])
set(gca, 'YTick', 0:0.5:1.5);
yl = ylim;
plot([1 1], yl, 'k--');
hp1 = fill([barX(1) barX(1) barX(end) barX(end)], [yl(2)*0.95 yl(2) yl(2) yl(2)*0.95], [1 1 1]*0.5, 'FaceAlpha', 0.5, 'EdgeColor', 'none');
hp1.FaceAlpha = 0.5;
uistack(hp1, 'bottom');
xlabel('time (s)');
ylabel('d''');
xlim(xRange);
offsetAxes();
spaceOutAxes();

% d' spatial map
p(1, 4).select();

cmap = parula(256);
x = imagesc(movData);
x.AlphaData = ~isnan(movData);
setImageAxis();
colormap(cmap);
axp = get(gca, 'Position');
c = colorbar;
c.Label.String = 'd''';

caxis([0 1.2]);

hold on;
for it3 = 1:length(contourListTransform)
  plot(contourListTransform{it3}(:,1),contourListTransform{it3}(:,2),'w');
end
for it3 = 1:length(contourListTransform)
  if(labelCenterListTransform{it3}(1) < 67 & labelCenterListTransform{it3}(2) < 67 & labelCenterListTransform{it3}(1) > 0 & labelCenterListTransform{it3}(2) > 0)
    text(labelCenterListTransform{it3}(1), labelCenterListTransform{it3}(2), areasNames{it3}, 'HorizontalAlignment', 'center', 'Color', [0 0 0]);
  end
end
axis off;
xlim([1 67]);
ylim([1 67]);  
set(gca, 'Position', axp-[0.025 0 0 0]);

set(findall(hFig,'-property','FontName'),'FontName', 'Arial');
set(findall(hFig,'-property','FontSize'),'FontSize', 8);