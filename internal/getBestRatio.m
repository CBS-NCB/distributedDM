function [rows, cols] = getBestRatio(numSquares, targetRatio)
  if(nargin < 2)
    targetRatio = 16/9;
  end
  squareList = 1:numSquares;
  prevMaxN = 0;
  bestRatioList = zeros(max(squareList), 3);
  for it0 = squareList
    if(prevMaxN >= it0)
      bestRatioList(it0, :) = bestRatioList(it0-1, :);
      continue;
    end
    ratioList = zeros(it0, 3);
    for itc = 1:it0
      nr = ceil(it0/itc);
      ratioList(itc, :) = [nr itc itc/nr];
    end
    [~, bestRatio] = min(abs(ratioList(:, 3)-targetRatio));
    bestRatioList(it0, :) = ratioList(bestRatio, :);
    prevMaxN = prod(ratioList(bestRatio, 1:2));
  end
  rows = bestRatioList(end, 1);
  cols = bestRatioList(end, 2);
end