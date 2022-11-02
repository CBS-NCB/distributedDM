function setImageAxis(varargin)
  if(nargin < 1)
    ax = gca;
  else
    ax = varargin{1};
  end
  axis(ax, 'tight', 'ij');
  set(ax, 'XTick', []);
  set(ax, 'YTick', []);
end