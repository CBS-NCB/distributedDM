clear all;
load('fig1d.mat');

hFig = createCenteredFigure('width',5,'height', 4);
hold on;
cmap = lines(1);
for it1 = 1:size(perfRI, 2)
  h = plot(binnedAnglesCenter, perfRI(:, it1)*100, '-', 'Color', [cmap(1, :) 0.5]);
  h.LineWidth = 0.5;
end
h = errorbar(binnedAnglesCenter, perfRM*100, perfRE*100, '-o', 'MarkerSize', 4, 'DisplayName', 'global', 'Color', cmap(1, :));
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


hFig = createCenteredFigure('width',5,'height', 4);
hold on;
cmap = lines(1);
for it1 = 1:size(perfTI, 2)
  h = plot(binnedAnglesCenter, perfTI(:, it1)*100, '-', 'Color', [cmap(1, :) 0.5]);
  h.LineWidth = 0.5;
end
h = errorbar(binnedAnglesCenter, perfTM*100, perfTE*100, '-o', 'MarkerSize', 4, 'DisplayName', 'global', 'Color', cmap(1, :));
h.MarkerFaceColor =  'k';

xlabel('Angle difference \langle|\theta_L|-|\theta_R|\rangle');
ylabel('Ignore (%)');
xlim([-90 90]);
ylim([0 20]);
set(gca, 'XTick', linspace(-90, 90, 5));
set(gca, 'YTick', linspace(0, 20,5));
spaceOutAxes(gca);
offsetAxes(gca, 50);
set(findall(hFig,'-property','FontName'),'FontName', 'Arial');
set(findall(hFig,'-property','FontSize'),'FontSize', 8);
