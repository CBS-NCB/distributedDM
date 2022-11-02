clear all;
close all;
clc;

pause(1);

%%

clear all;
% close all;
clc;

load ./M150617/AnalyzedRet.mat

% h = imshow(imresize(abs(CmplxMaps{1}{1}),0.5,'nearest'),[]);
% [I2 rect] = imcrop(h);
% save cropRegion.mat rect;

load cropRegion.mat;

CmplxMaps{1}{1} = imresize(CmplxMaps{1}{1},0.5,'nearest');
CmplxMaps{1}{2} = imresize(CmplxMaps{1}{2},0.5,'nearest');
CmplxMaps{2}{1} = imresize(CmplxMaps{2}{1},0.5,'nearest');
CmplxMaps{2}{2} = imresize(CmplxMaps{2}{2},0.5,'nearest');

P = round(rect);

CmplxMaps{1}{1} = CmplxMaps{1}{1}(P(2):P(2)+P(4),P(1):P(1)+P(3));
CmplxMaps{1}{2} = CmplxMaps{1}{2}(P(2):P(2)+P(4),P(1):P(1)+P(3));
CmplxMaps{2}{1} = CmplxMaps{2}{1}(P(2):P(2)+P(4),P(1):P(1)+P(3));
CmplxMaps{2}{2} = CmplxMaps{2}{2}(P(2):P(2)+P(4),P(1):P(1)+P(3));

ampThresh = 0;
pixpermm = 157.8/2;
datafeatures = func_ExtractFeatures(CmplxMaps,visual_field,ampThresh,pixpermm);


%%
addpath /media/nev/Data/DRTOOLBOX;
addpath /media/nev/Data/DRTOOLBOX/gui;
addpath /media/nev/Data/DRTOOLBOX/techniques;

nodims = size(datafeatures,2);
knearest = 12;
mappedX = compute_mapping(datafeatures,'Isomap', nodims,'k',knearest);  
matFileName = sprintf('./ISOMAPresults/isomapTensorAll_%dDim_k%d.mat',nodims,knearest);
save(matFileName,'mappedX');





