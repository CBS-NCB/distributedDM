function ax = offsetAxes(ax, varargin)
% thanks to Pierre Morel, undocumented Matlab
% and https://stackoverflow.com/questions/38255048/separating-axes-from-plot-area-in-matlab
%
% by Anne Urai, 2016
% Modified by Javier Orlandi, 2019

if(~exist('ax', 'var'))
  ax = gca;
end
ax.XAxis.Color = 'k';
ax.YAxis.Color = 'k';
set(ax, 'TickDir', 'out');
box(ax, 'off');
set(ax.Parent,'Color','w');

if(nargin > 1)
  factor = varargin{1};
else
  factor = 20;
end
  % modify the x and y limits to below the data (by a small amount)
  XT = ax.XTick;
  YT = ax.YTick;
  ax.XTickMode = 'manual';
  ax.YTickMode = 'manual';

  

  ax.XLim(1) = ax.XLim(1)-(ax.XLim(2)-ax.XLim(1))/factor;
  ax.YLim(1) = ax.YLim(1)-(ax.YLim(2)-ax.YLim(1))/factor;
  %ax.XLim(1) = ax.XLim(1)-(ax.XTick(2)-ax.XTick(1))/factor;
  %ax.YLim(1) = ax.YLim(1)-(ax.YTick(2)-ax.YTick(1))/factor;
  % this will keep the changes constant even when resizing axes
  addlistener(ax, 'MarkedClean', @(obj,event)resetVertex(ax, factor));
end

function resetVertex (ax, fac)
  factor = 1/fac;
  %ax.XTickMode = 'manual';
  % extract the x axis vertext data
  % X, Y and Z row of the start and end of the individual axle.
  %ax.XRuler.Axle.VertexData
  %ax.XRuler.Axle.VertexData(1,1) = min(get(ax, 'Xtick'));
  ax.XRuler.Axle.VertexData(1,1) = (ax.XLim(1)+ax.XLim(2)*factor)/(1+factor);
  % repeat for Y (set 2nd row)
  %ax.YRuler.Axle.VertexData(2,1) = min(get(ax, 'Ytick'));
  ax.YRuler.Axle.VertexData(2,1) = (ax.YLim(1)+ax.YLim(2)*factor)/(1+factor);
  %ax.XRuler.Axle.VertexData
end