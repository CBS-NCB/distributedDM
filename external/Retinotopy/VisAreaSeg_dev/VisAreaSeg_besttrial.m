% Visual Area Segmentation (Garret 2014 J Neurosci) MA 20150507
clearvars; close all; clc; clear global

%% Declarations and setting required directories
global DIRS VDAQ
SetBaseDir;
% addpath(fullfile(DIRS.baseDir,'\LAB_CODE\MATLAB\MOUSE\WIDEFIELD\Main_WideField'));
indicator = 'GCAMP'; animal= 'M151120_MA'; iseries = 2; Expts = [1 2]; Serv = 2; 
func_DefaultDirs(indicator, animal, iseries, Serv);
CameraInfo.HardwareBinning  = 1; % binning horz and vert
CameraInfo.MmPerCameraPix   = 0.0063 * CameraInfo.HardwareBinning; 

%% Load the SNR analysis and label the trials to actually import
load('C:\Users\Mike\Documents\201512_widefieldQC\res_20151211_2.mat'); % 'X' is loaded
% match the day in question
snrsess = {X.animal}; matchsess = ismember(snrsess,animal);
snrseries = [X.iseries]; matchseries = ismember(snrseries,iseries);
VDAQ.X = X(matchsess&matchseries);

%% here, Load the data and do basic prep using Main_WideField.m for all the iexps
for idx  = 1:length(Expts)
    iexp  = Expts(1,idx);
    load([DIRS.ScreenInfoDir '_' num2str(iexp)]); clc;   % Load ScreenInfo
    VDAQ.Xidx = idx; % if "Expts" doesn't start on 1, we need a way to know
    
    RepFlag = 0; % to keep single trial data (not just averaged)
    resFac = .8; 
    Main_WideField_besttrial; % Run the initial analyses
    
    [AbsMaps{idx}, AngleMaps{idx}, CmplxMaps{idx}] = TensorFrequency([], myfreqs ); % Get Azm and elv maps
    VDAQ_ALL(idx)= VDAQ; % Keep both VDAQs in workspace in case you need it.
    
    % Calculate visual field (screen x and y in degrees of visual angle)
    monitor_distance = myScreenInfo.Dist;
    if strcmp(p.xfile,'stimMarshel.x'),
        if  p.pars(13,1)==1 % ori 90 or 0 (if 1 or 2)
            screen_fractionX = (p.pars(3,1)-p.pars(2,1))/360;
            L = myScreenInfo.Xmax*myScreenInfo.PixelSize*screen_fractionX; % in cm
        elseif p.pars(13,1)==2,
            screen_fractionY = (p.pars(3,1)-p.pars(2,1))/360;
            L = myScreenInfo.Ymax*myScreenInfo.PixelSize*screen_fractionY; % in cm
        end    
    elseif strcmp(p.xfile,'stimKalatsky.x'),
        if  p.pars(end,1)==1 % ori 90 or 0 (if 1 or 2)
            screen_fractionX = (p.pars(3,1)-p.pars(2,1))/360;
            L = myScreenInfo.Xmax*myScreenInfo.PixelSize*screen_fractionX; % in cm
        elseif p.pars(end,1)==2,
            screen_fractionY = (p.pars(3,1)-p.pars(2,1))/360;
            L = myScreenInfo.Ymax*myScreenInfo.PixelSize*screen_fractionY; % in cm
        end
    end
    visual_field(idx) = 2*atan(L/(2*monitor_distance))*180/pi;% in visual angle
end

% Apply eye position (meridians) with respect to the screen
EyeY = 0.; 
EyeX = 0;

%SingRepVFS_VDAQALL


%% Choose ROIs if needed (load the analyzedret if needed) and get Direct and Reverse maps for altitude and azimuth
% load([DIRS.SaveDir 'AnalyzedRet'])
% [Roi,~] = ClipROI(CmplxMaps{1});
% for i =1:2
%     for j =1:2
%         CmplxMapsRoi{i}{j} = CmplxMaps{i}{j}(Roi(2):Roi(2)+Roi(4), Roi(1):Roi(1)+Roi(3));
%     end
% end
% ang0 = angle(CmplxMapsRoi{1}{1});
% ang2 = angle(CmplxMapsRoi{1}{2});
% ang1 = angle(CmplxMapsRoi{2}{1});
% ang3 = angle(CmplxMapsRoi{2}{2});

%% Direct and Reverse maps for altitude and azimuth
ang0 = angle(CmplxMaps{1}{1});
ang2 = angle(CmplxMaps{1}{2});
ang1 = angle(CmplxMaps{2}{1});
ang3 = angle(CmplxMaps{2}{2});

%% Segmentation (Ian Nauhaus)
%Find delay as the angle between the vectors
delay_hor  = angle( exp(1i*ang0)+exp(1i*ang2) );
% Other corrections on the delaymaps
% angind = delay_hor > 0;
% delay_hor(angind) = delay_hor(angind)-2*pi;
% subplot(223); imagesc(delay_hor); axis image; caxis([0 pi])

delay_vert = angle( exp(1i*ang1)+exp(1i*ang3) );

%The delay always pushes the vectors counter clockwise. (This is simply mod(val, pi), MA20150727)
delay_hor  = delay_hor  + (pi/2)*(1-sign(delay_hor));
delay_vert = delay_vert + (pi/2)*(1-sign(delay_vert));

%Use delay vector to calculate retinotopy.
map_hor  = 0.5*(angle(exp(1i*(ang0-delay_hor)))  - angle(exp(1i*(ang2-delay_hor))));
map_vert = 0.5*(angle(exp(1i*(ang1-delay_vert))) - angle(exp(1i*(ang3-delay_vert))));

map_hor  = (map_hor  /(2*pi)) * visual_field(1);
map_vert = (map_vert /(2*pi)) * visual_field(2);

%% Smooth the maps if needed and then get the visual field sign map
% GF = fspecial('gaussian',size(map_hor), 3);
% GF = GF/sum(GF(:));
% map_hor = ifft2( fft2(map_hor).*abs(fft2(GF)) );
% map_vert = ifft2( fft2(map_vert).*abs(fft2(GF)) );

[VFS, VFS_thr] = getVisSign(map_hor, map_vert);

%% Save the file for segmentation analysis
% save([DIRS.SaveDir 'AnalyzedRet'],'CmplxMaps','EyeX', 'EyeY','visual_field')








%% Check the delay maps and remap if needed or if delay maps are strange
figure; 
subplot(221); imagesc(delay_hor); axis image; caxis([0 pi])
subplot(222); imagesc(delay_vert); axis image; caxis([0 pi])

delay_hor_1 = delay_hor;
angind = delay_hor_1 < pi/1.5;
delay_hor_1(angind) = abs(delay_hor_1(angind) - pi);
subplot(223); imagesc(delay_hor_1); axis image; caxis([0 pi])

delay_vert_1 = delay_vert;
angind = delay_vert_1 < pi/1.5;
delay_vert_1(angind) = abs(delay_vert_1(angind) - pi);
subplot(224); imagesc(delay_vert_1); axis image; caxis([0 pi])

%Use remapped delays vector to calculate retinotopy.
map_hor_1  = 0.5*(angle(exp(1i*(ang0-delay_hor_1)))  - angle(exp(1i*(ang2-delay_hor_1))));
map_vert_1 = 0.5*(angle(exp(1i*(ang1-delay_vert_1))) - angle(exp(1i*(ang3-delay_vert_1))));

map_hor_1  = (map_hor_1  /(2*pi)) * visual_field(1);
map_vert_1 = (map_vert_1 /(2*pi)) * visual_field(2);

[VFS1, VFS_thr1] = getVisSign(map_hor_1, map_vert_1);

%% variances of each repeat
% for istim=1:p.nstim
% %     MyMin = prctile( VDAQ_ALL(1).repeatVar{istim}(:), 1);
% %     MyMax = prctile( VDAQ_ALL(1).repeatVar{istim}(:), 99);
% %     figure; 
% %     for irep=1:p.nrepeats, 
% %         subplot(4,5,irep); imshow(fliplr(rot90(VDAQ_ALL(1).repeatVar{istim}(:,:,irep),3)),[MyMin MyMax]); axis image
% %     end
% %     VarofVars = (var(VDAQ_ALL(1).repeatVar{istim},[],3));
% %     figure; imshow(fliplr(rot90(VarofVars,3)),[]); axis image
% 
%     MyMin = prctile( VDAQ_ALL(2).repeatVar{istim}(:), 1);
%     MyMax = prctile( VDAQ_ALL(2).repeatVar{istim}(:), 99);
%     figure; 
%     for irep=1:20, 
%        subplot(4,5,irep); imshow(fliplr(rot90(VDAQ_ALL(2).repeatVar{istim}(:,:,irep),3)),[MyMin MyMax]); axis image
%     end
%     %VarofVars = (var(VDAQ_ALL(2).repeatVar{1},[],3));
%     %figure; imshow(fliplr(rot90(VarofVars,3)),[]); axis image
% end
% 
% %% For Dmitry
% % superimpose VFS and Segmeted Areas 
% AreasOrigSz = imresize(VisAreas,size(ang0),'nearest');
% BW = bwlabel(AreasOrigSz,4);
% figure
% imshow(fliplr(rot90(VDAQ_ALL(2).repeatVar{istim}(:,:,irep),3)),[MyMin MyMax])
% hold on;
% [C,H]= imcontour(BW,[1 1],'r');
% set(H,'linewidth',2)

%save([DIRS.SaveDir 'VisAreas_mouse39.mat'],'AreasOrigSz','BW');






