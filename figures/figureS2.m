%%
% CD to this folder
addpath(genpath(fullfile(pwd,'helperFunctions')));

%%
clearvars;
load('data/figs2.mat');

hFig = createCenteredFigure('width', 20, 'height', 20);
p = panel();
p.margin = 1;
p.pack(10, 10);
regionNames = {'V1', 'PM', 'AM', 'A', 'SSt', 'RL', 'SSb', 'AL', 'VL', 'RS'};

for it1 = 1:10

  for it2 = 1:length(fullMat{it1})
    p(it1,it2+1).select();
    hold off;

    imagesc(fullMat{it1}{it2});
    if(it2 == 1)
      ylabel(regionNames{it1});
      ax = gca;
      ax.YLabel.Rotation = 0;
      ax.YLabel.HorizontalAlignment = 'right';
    end
    caxis([0.05, 0.95]);
    hold on;
    for it3 = 1:length(contourListTransform)
      plot(contourListTransform{it3}(:,1),contourListTransform{it3}(:,2),'w');
    end
    setImageAxis();
    xlim([1 67]);
    ylim([1 67]);
    colormap(parula(256));
    %title(regionNames{curComp});  
  end
end
set(findall(hFig,'-property','FontName'),'FontName', 'Arial');
set(findall(hFig,'-property','FontSize'),'FontSize', 7);

