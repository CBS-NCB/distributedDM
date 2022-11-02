function [xRange, yRange] = spaceOutAxes(varargin)
% SPACEOUTAXES Adds a 2.5% to the axis limits (for breathing purposes)
%
% USAGE:
%   spaceOutAxes();
%
% INPUT arguments:
%   ax - axes of choice (if empty will use the current axes)
%
% OUTPUT arguments:
%   xRange, yRange - the new x and  y limits
%
% EXAMPLE:
%   [xRange, yRange] = spaceOutAxes();
%
% Copyright (C) 2019, Javier G. Orlandi <javierorlandi@javierorlandi.com>,
% Benuccilab (http://benuccilab.brain.riken.jp)

  if(nargin == 0)
    ax = gca;
  else
    ax = varargin{1};
  end
  xRange = xlim(ax);
  xlim(ax, xRange+[-1 1]*diff(xRange)*0.025);
  yRange = ylim(ax);
  ylim(ax, yRange+[-1 1]*diff(yRange)*0.025);
end