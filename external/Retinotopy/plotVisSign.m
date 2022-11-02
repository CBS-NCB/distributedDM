function varargout = plotVisSign(VSgnMap, pixpermm)
% INPUTS:
%     map_horz - Map of horizontal retinotopy
%     map_vert - Map of vertical retinotopy
%     pixpermm = 1mm / pixSize(in mm) of the retinotopy images
% OUTPUTS:
%     VisSgnMap is a the visual field sign map
%     VS_Thr is threshold of VisSgnMap 

%%
if nargin < 2
    pixpermm = 1/0.0063; %The size of each pixel is 0.0063 in PCO images
end

%% Plot absolute horz and vert retinotopic maps
mmperpix = 1/pixpermm;
xdom = (0:size(VSgnMap,2)-1)*mmperpix;
ydom = (0:size(VSgnMap,1)-1)*mmperpix;

%%
AbsVFS  = abs(VSgnMap);
Thresh  = 1.5 * std(VSgnMap(:));
VS_Thr  = (sign(AbsVFS - Thresh/2) + 1)/2;  %threshold VSgnMap at (+-) 1.5 std

id      = find(VS_Thr);
imdum   = VS_Thr.*VSgnMap; 
imdum(id)= imdum(id)+1.1;


%% Plotting thresholded patches of VSgnMap
figure;
ploteccmap(imdum, [.1 2.1], pixpermm);
colorbar % colorbar off
axis image
StrThr = num2str(Thresh);
title(['Threshold: ' StrThr(1:4)])