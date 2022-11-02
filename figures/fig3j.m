clearvars;
load('fig3j.mat');

%% Final figure

hFig = createCenteredFigure('width', 5, 'height', 5);
cmap = cblindmap();
cmap = cmap(1:end,:);
hold on;
  
%%% The projection
for it0 = 1:4
  switch it0
    case 1
      attMult = 1;
      diffMult = 1;
    case 2
      attMult = 1;
      diffMult = 0;
    case 3
      attMult = 0;
      diffMult = 1;
    case 4
      attMult = 0;
      diffMult = 0;
  end
  opacity = (1+diffMult)/2;
  h1 = plot3(x(:, (it0-1)*2+1), y(:, (it0-1)*2+1), z(:, (it0-1)*2+1), '-', 'Color', cmap(1, :)*opacity);
  h2 = plot3(x(:, (it0-1)*2+2), y(:, (it0-1)*2+2), z(:, (it0-1)*2+2), '-', 'Color', cmap(1+1, :)*opacity);
end

yl = ylim;
plot([0 0], ylim, 'k--');
set(gca, 'XTick', 0:2);
xlabel('time (s)');
ylabel('Projection (a.u.)');
xlim([-0.5 0.5]);
set(gca,'XTick',-0.5:0.5:1);
set(gca,'YTick',-2:1:2);
ylim([-2 2]);
zlim([0 1.2]);
grid on;
set(gca,'ZTick',0:0.25:1);
view([34.5 30.5]);
