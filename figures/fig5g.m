clearvars;
load('fig5g.mat');

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
