function baseImg = getSVDcomponent(dataSVD, compIndex)

  baseImg = nan(size(dataSVD.globalMask));
  validPixels = dataSVD.valid;
  baseImg(validPixels) = 0;
  validPixels = find(baseImg' == 0);

  baseImg = nan(size(dataSVD.globalMask))';
  if(compIndex == 0)
    baseImg(validPixels) = dataSVD.means;
  else
    baseImg(validPixels) = dataSVD.coeff(:, compIndex);
  end
  baseImg = baseImg';
end
