%%
clearvars;
load('figs1b.mat');

%%
hFig = createCenteredFigure('width', 8, 'height', 6);
p = panel();
p.margin = 2;
p.pack(2, 5);
regionNames = {'V1', 'PM', 'AM', 'A', 'SSt', 'RL', 'SSb', 'AL', 'VL', 'RS'};
curComp = 0;

for it1 = 1:2
  for it2 = 1:5
    curComp = curComp + 1;
    p(it1,it2).select();
    hold off;
    imagesc(1-baseImgFull(:, :, curComp));
    %imagesc(baseImgFull(:, :, curComp));
    hold on;
    axis equal;
    for it3 = 1:length(contourListTransform)
      plot(contourListTransform{it3}(:,1),contourListTransform{it3}(:,2),'k');
    end
    setImageAxis();
    xlim([1 67]);
    ylim([1 67]);
    colormap([1 1 1;parula]);
    title(regionNames{curComp});  
  end
end
set(findall(hFig,'-property','FontName'),'FontName', 'Arial');
set(findall(hFig,'-property','FontSize'),'FontSize', 7);

