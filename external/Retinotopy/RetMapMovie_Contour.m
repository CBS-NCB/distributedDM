%% need to get the patches during debugging mode (should get these info as output of the VFS toolbox)
clear all; close all; clc;

%% Declarations
global DIRS VDAQ
SetBaseDir;
% addpath(fullfile(DIRS.baseDir,'\LAB_CODE\MATLAB\MOUSE\WIDEFIELD\Main_WideField'));
indicator = 'INTR'; animal= 'M150527_MA_TT'; iseries =  1; iexp = 1; 
func_DefaultDirs(indicator, animal, iseries);

load([DIRS.SaveDir, 'AnalyzedRet.mat']);

%% Load Tensor
resFac = .8; ResizeF = num2str(resFac);
if resFac<1, cResizeF = sprintf('0%1.0f',round(resFac*10)); end
Main_WideField;

%%
% se1 = strel('disk',5);
% roibw = imdilate(roibw,se1);
% roibw = imdilate(roibw,se1);

%%
avrAmp = zeros(size(CmplxMaps{1}{1}));
for j=1:2
    for i=1:2
        temp = mat2gray(log(abs(CmplxMaps{1,j}{i,1})+1));
        avrAmp = avrAmp + temp/4;
    end
end
roibw = double(avrAmp>0.15);
%%
for j=1:2
    for i=1:2
        tempAmp = abs(CmplxMaps{1,j}{i,1});
        tempPhi = angle((CmplxMaps{1,j}{i,1}));
        Complex_scaledwith_Amp{1,j}{i,1} = avrAmp.* exp( sqrt(-1)*double(roibw).*tempPhi );
    end
end
bpViewComplexMaps(Complex_scaledwith_Amp{1,2})

%% Converting complex maps to degree of visual angle
Divs_Azm_tmp = Complex_scaledwith_Amp{1,1}{1,1} ./ Complex_scaledwith_Amp{1,1}{2,1};  
Divs_Elv_tmp = Complex_scaledwith_Amp{1,2}{1,1} ./ Complex_scaledwith_Amp{1,2}{2,1}; 

Azm_tmp = angle(Divs_Azm_tmp) /(2*pi)*visual_field(1,1);
Elv_tmp = angle(Divs_Elv_tmp) /(2*pi)*visual_field(1,2);

%% Get the visual field sign map (VFS) and its threshold (VFS_thr)
[VFS_tmp, VFS_thr_tmp] = getVisSign(Azm_tmp - EyeX, Elv_tmp - EyeY);
VFS_thr_tmp_labeled = bwlabel(VFS_thr_tmp);
for lbind = 1:max(VFS_thr_tmp_labeled(:))
    lbroi = VFS_thr_tmp_labeled == lbind;
    
    lbroiprop = regionprops(lbroi,'area');
    roiAreas(lbind,:) = [lbind lbroiprop.Area*mean(avrAmp(lbroi))];
end
[roiAreasSorted,roiAreasInd] = sortrows(roiAreas,2);
removeIndexes = roiAreasSorted(1:end-10,1);

for rmind = 1:length(removeIndexes)
    z = VFS_thr_tmp_labeled == removeIndexes(rmind);
    VFS_thr_tmp_labeled(z) = 0;    
end
VFS_thr_tmp_new = double(VFS_thr_tmp_labeled > 0);
VFS_tmp_new = VFS_tmp.*VFS_thr_tmp_new;

plotVisSign(VFS_tmp_new);

%% Superimpose the area contours on the VDAQ.tensor movies 
figure
for stimNum =1:2
    mymin   = prctile(VDAQ.tensor{stimNum}(:), 2.5);
    mymax   = prctile(VDAQ.tensor{stimNum}(:), 97.5);
    for frm=1:size(VDAQ.tensor{1},3)
        drawnow
        subplot(111)
        imagesc(squeeze(VDAQ.tensor{stimNum}(:,:,frm)), [mymin, mymax])
%         hold on; [CCC,HHH] = imcontour(VFS_thr_tmp_new, 1,'r-'); set(HHH,'LineWidth',2)
%         title(['Time: ' num2str(frm*(1/floor(VDAQ.FrameRate))) ' Sec'], 'fontsize', 15)
        colormap('gray'); colorbar;
        axis image
        TheMovie(stimNum, frm) = getframe;
    end
    filename    = [animal '_' num2str(iseries) '_' num2str(iexp) '_' num2str(stimNum) '_Contours.avi'];
    movie2avi(TheMovie(stimNum,:), [DIRS.SaveDir filename], 'compression','None', 'fps', floor(VDAQ.FrameRate));
end

%% 












 
