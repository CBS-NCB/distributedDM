clearvars;
load('fig4e.mat');

%% Final figure

dataName = {'Dorsal stream', 'Posterior stream', 'Ventral stream', 'Somatosensory', 'Retrosplenial'};
labels = {'Stimulus', 'Wheel mov', 'Saccades', 'Attention', 'Choice pre', 'Choice post'};
pairList = [1 2;1 3; 2 3; 4 5];

hFig = createCenteredFigure('width', 21, 'height', 9);
p = panel();
p.margin = 12;
p.pack([0.75 0.25], [[1, 1, 1, 1]*0.235, 0.004]);
cmap = cblindmap;
cmap = cmap([1 2 3 6 5 4 7], :);
ax = [];
for it0 = 1:4
  p(1, it0).select();
  hold on;
  h = [];
  pair_x = pairList(it0, 1);
  pair_y = pairList(it0, 2);
  x = data{pair_x};
  y = data{pair_y};
  for it1 = 1:6
    sub_x = x(:, it1);
    sub_y = y(:, it1);
    %curData = squeeze(fullData(:, :, it1, itSel));
    plot(sub_x, sub_y, '.', 'Color', cmap(it1, :), 'MarkerFace', cmap(it1, :));
    h = [h; errorbar(mean(sub_x), mean(sub_y), -std(sub_y)/sqrt(7), std(sub_y)/sqrt(7), -std(sub_x)/sqrt(7), std(sub_x)/sqrt(7), 'o', 'Color', cmap(it1, :), 'MarkerSize', 4, 'MarkerFace', cmap(it1, :))];
  end
  uistack(h, 'top');
  xlabel(dataName{pair_x});
  ylabel(dataName{pair_y});
  axis equal square tight;
  xlim([0 1.5]);
  ylim([0 1.5]);
  hold on;
  plot([0 1.5], [0 1.5], 'k--');
  offsetAxes();
  spaceOutAxes();
  ax = [ax; gca];
end
l = legend(h, labels);
l.Box = 'off';
l.Position(1) = l.Position(1) +0.13;