function [runningVar, runningN] = onlineVariance(x)
% ONLINEVARIANCE Online variance computation
%
% USAGE:
%   runningVar = onlineVariance(x);
%
% INPUT arguments:
%   x - new data point (empty to initialize a new run)
%
% OUTPUT arguments:
%   runningVar - current variance value
%   n - number of samples
%
% EXAMPLE:
%   onlineVariance([])
%   x = zeros(1000, 1);
%   n = zeros(1000, 1);
%   for it = 1:1000
%     [x(it) n(it)] = onlineVariance(randn);
%   end
%   figure;
%   plot(n, x);
%   xlabel('iteration');
%   ylabel('Running variance');

persistent M S n;
if(isempty(x) || isempty(M) || isempty(S) || isempty(n))
  M = 0;
  S = 0;
  n = 0;
else
  n = n + 1;
  oldM = M;
  M = M + (x-M)/n;
  S = S + (x-M).*(x-oldM);
end
runningVar = S/(n-1);
runningN = n;

