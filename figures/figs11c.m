%%
clearvars;
load('figs11c.mat');

%%

regionNames = {'V1', 'PM', 'AM', 'A', 'SSt', 'RL', 'SSb', 'AL', 'L', 'RS'};
hFig = createCenteredFigure('width', 21, 'height', 10);
p = panel();
p.margin = 10;
p.pack(1, [0.3, 0.3 0.4]);
p(1,1).select();

h = imagesc(imgPre);
set(h, 'AlphaData', ~isnan(imgPre));
colorbar; 
title('before movement');
tf = [];
hold on;
for it3 = 1:length(contourListTransform)
  plot(contourListTransform{it3}(:,1),contourListTransform{it3}(:,2),'w');
end
for it3 = 1:length(contourListTransform)
  if(labelCenterListTransform{it3}(1) < 67 & labelCenterListTransform{it3}(2) < 67 & labelCenterListTransform{it3}(1) > 0 & labelCenterListTransform{it3}(2) > 0)
    tf =[tf; text(labelCenterListTransform{it3}(1), labelCenterListTransform{it3}(2), regionNames{it3}, 'HorizontalAlignment', 'center', 'Color', [0 0 0])];
  end
end
tf = [];
setImageAxis();
axis equal;
axp = get(gca, 'Position');
c = colorbar;
c.Label.String = 'd''';
axis off;
xlim([1 67]);
ylim([1 67]);  
set(gca, 'Position', axp-[0.025 0 0 0]);


p(1,2).select();
h = imagesc(imgPost);
set(h, 'AlphaData', ~isnan(imgPost));
colorbar;
title('after movement');
hold on;
for it3 = 1:length(contourListTransform)
  plot(contourListTransform{it3}(:,1),contourListTransform{it3}(:,2),'w');
end
for it3 = 1:length(contourListTransform)
  if(labelCenterListTransform{it3}(1) < 67 & labelCenterListTransform{it3}(2) < 67 & labelCenterListTransform{it3}(1) > 0 & labelCenterListTransform{it3}(2) > 0)
    tf =[tf; text(labelCenterListTransform{it3}(1), labelCenterListTransform{it3}(2), regionNames{it3}, 'HorizontalAlignment', 'center', 'Color', [0 0 0])];
  end
end
setImageAxis();
axis equal;
axp = get(gca, 'Position');
c = colorbar;
c.Label.String = 'd''';
axis off;
xlim([1 67]);
ylim([1 67]);  
set(gca, 'Position', axp-[0.025 0 0 0]);


set(findall(hFig,'-property','FontName'),'FontName', 'Arial');
set(findall(hFig,'-property','FontSize'),'FontSize', 8);
