clearvars;
load('fig2a.mat');

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