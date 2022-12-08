%%
clearvars;
load('figs1a.mat');

%% The same with the additional transformation to the 67 67
hFig = createCenteredFigure('width', 5, 'height', 5);
hold on;
for it3 = 1:length(contourListTransform)
  plot(contourListTransform{it3}(:,1),contourListTransform{it3}(:,2),'k');
  if(labelCenterListTransform{it3}(1) < 67 & labelCenterListTransform{it3}(2) < 67 & labelCenterListTransform{it3}(1) > 0 & labelCenterListTransform{it3}(2) > 0)
    text(labelCenterListTransform{it3}(1), labelCenterListTransform{it3}(2), areasNamesOrig{it3}, 'HorizontalAlignment', 'center');
  end
end
setImageAxis();
xlim([1 67]);
ylim([1 67]);
box on;

for it1 = 1:7
  ax = gca;
  ax.ColorOrderIndex = 1;
  for it2 = 1:3
    hb = plot(area_x{it1}{it2}, area_y{it1}{it2}, 'DisplayName', fnames{it2});
  end
end

set(findall(hFig,'-property','FontName'),'FontName', 'Arial');
set(findall(hFig,'-property','FontSize'),'FontSize', 7);

