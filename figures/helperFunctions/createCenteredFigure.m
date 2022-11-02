function hFig = createCenteredFigure(varargin)
% CREATECENTEREDFIGURE Creates a new figure centered on the screen
%
% USAGE:
%   hFig = createCenteredFigure(varargin)
%
% INPUT optional arguments ('key' followed by its value):
%   width - figure width (default 29.7 cm)
%   height - figure height (default 21 cm)
%   units - units for the figure (including width and height) Default (centimeters)
%   visible - If the new figure should be visible (default true)
%   useMultipleMonitors - If the figure should be centered across monitors (default false, If false, it will use the current monitor).
%   maximized - Will maximize the figure after drawing it
%   parent - If the figure should be centered across monitors (default false, If false, it will use the current monitor).
%   ppi - PPI of your monitor (default: ScreenPixelsPerInch) - should be 96 on windows
%
% OUTPUT arguments:
%   hFig - figure handle
%
% EXAMPLE:
%   hFig = createCenteredFigure('width', 10, 'height', 20)
%
% Copyright (C) 2016-2018, Javier G. Orlandi <javiergorlandi@gmail.com>
 
params.units = 'centimeters';
params.width = 29.7; % A4
params.height = 21.05;
params.ppi = []; % PPI adjustment (in case your monitor is not 96 ppi)
params.visible = true;
params.useMultipleMonitors = false;
params.maximized = false; % Will maximize at the end
% Parse them
params = parse_pv_pairs(params, varargin);

% Let's create a dummy figure to get conversion ratios
h = figure('Visible', 'off');
h.Position(3:4) = [1000 1000];
h.Units = params.units;
ratio = 1000./h.Position(3:4); % In px/units (2 ratios in case some weird pixel to cm stuff)
if(~isempty(params.ppi))
  ratio = params.ppi/get(groot,'ScreenPixelsPerInch')*ratio;
end
delete(h);
widthPixels = ratio(1)*params.width;
heightPixels = ratio(2)*params.height;

if(ismac)
  BORDERSIZE = 150;
else
  BORDERSIZE = 50;
end
% Get absolute available area
monPos = get(0,'MonitorPositions');
currentMonitor = 1;
if(params.useMultipleMonitors)
  monMinX = min(monPos(:,1));
  monMaxX = max(monPos(:,1)+monPos(:,3)-1);
  monMinY = min(monPos(:,2));
  monMaxY = max(monPos(:,2)+monPos(:,4)-1);  
else

  % Get current matlab monitor
  desktop = com.mathworks.mde.desk.MLDesktop.getInstance;
  desktopMainFrame = desktop.getMainFrame;
  % That +9 is so weird....
  mainX = desktopMainFrame.getLocation.x+9;
  mainY = desktopMainFrame.getLocation.y+9;
  if(mainX < 1)
    mainX = 1;
  end
  if(mainY < 1)
    mainY = 1;
  end
 
  currentMonitor = 1;
  if(size(monPos, 1) > 1)
    % If there's more than one monitor, check in which one we are
    for it = 1:(size(monPos, 1))
      if(mainX >= monPos(it, 1) && mainX < (monPos(it, 1)+monPos(it,3)) && mainY >= monPos(it, 2) && mainY < (monPos(it, 2)+monPos(it,4)))
        currentMonitor = it;
        break;
      end
    end
    monMinX = monPos(currentMonitor,1);
    monMaxX = monPos(currentMonitor,1)+monPos(currentMonitor,3)-1;
    monMinY = monPos(currentMonitor,2);
    monMaxY = monPos(currentMonitor,2)+monPos(currentMonitor,4)-1;
  else
    monMinX = min(monPos(:,1));
    monMaxX = max(monPos(:,1)+monPos(:,3)-1);
    monMinY = min(monPos(:,2));
    monMaxY = max(monPos(:,2)+monPos(:,4)-1);
  end
end

% Center it on the current monitor
pos = [monPos(currentMonitor,1)+round((monPos(currentMonitor,3)-widthPixels)/2), monPos(currentMonitor, 2)+round((monPos(currentMonitor,4)-heightPixels)/2), widthPixels, heightPixels];
% If it still doesn't fit, change width and height so it fits
if(pos(1)+pos(3) + BORDERSIZE > monMaxX || pos(2) + pos(4) + BORDERSIZE > monMaxY || pos(1) < monMinX || pos(2) < monMinY)
  if(widthPixels > monPos(1,3))
    widthPixels = monPos(1,3)-BORDERSIZE;
  end
  if(heightPixels > monPos(1,4))
    heightPixels = monPos(1,4)-BORDERSIZE;
  end
  fprintf('Could not fit figure on the screen. Resizing to: %dx%d px.\n', widthPixels, heightPixels);
  % If it still doesn't fit, you are pretty much screwed
  pos = [monPos(1,1)+round((monPos(1,3)-widthPixels)/2), monPos(1, 2)+round((monPos(1,4)-heightPixels)/2), widthPixels, heightPixels];
end

% The actual figure
hFig = figure('Visible', params.visible, 'units', 'pixels');
hFig.Color = 'w';
hFig.Position = pos;
if(params.maximized)
  hFig.WindowState = 'maximized';
end
% Do not. Better stick to the default wich is pixels
%hFig.Units = params.units; % Convert to appropiate units
