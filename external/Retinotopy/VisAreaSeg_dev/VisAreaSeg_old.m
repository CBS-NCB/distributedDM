% Visual Area Segmentation (Garret 2014 J Neurosci) MA 20150507
% clearvars; close all; clc; 
clear global

%% Declarations and setting required directories
global DIRS VDAQ
SetBaseDir;

indicator = 'GCAMP'; animal= 'M160727_MA'; iseries = 1; Expts = [1 2]; Serv = 2; 

func_DefaultDirs(indicator, animal, iseries, Serv);
CameraInfo.HardwareBinning  = 1; % binning horz and vert
CameraInfo.MmPerCameraPix   = 0.0063 * CameraInfo.HardwareBinning; 
multsess = false;

%% here, Load the data and do basic prep using Main_WideField.m for all the iexps
for idx  = 1:length(Expts)
    iexp  = Expts(1,idx);
    if exist('oldExpts','var'), oldiexp = oldExpts(1,idx); end;
    load([DIRS.ScreenInfoDir '_' num2str(iexp)]);    % Load ScreenInfo
    
    RepFlag = 0; % to keep single trial data (not just averaged)
    Main_WideField; % Run the initial analyses

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
EyeY = 0; 
EyeX = 0;

%SingRepVFS_VDAQALL

%% Direct and Reverse maps for altitude and azimuth
ang0 = angle(CmplxMaps{1}{1});
ang2 = angle(CmplxMaps{1}{2});
ang1 = angle(CmplxMaps{2}{1});
ang3 = angle(CmplxMaps{2}{2});

%% phase elev/azim/sign maps (orig author: Ian Nauhaus)
% find delay as the angle between the vectors
delay_hor  = angle( exp(1i*ang0)+exp(1i*ang2) );
delay_vert = angle( exp(1i*ang1)+exp(1i*ang3) );

% make delay go from 0 to pi and 0 to pi, instead of 0 to pi and 0 to -pi;
% delay can't be negative. if the delay vector is in the bottom two 
% quadrants, it is assumed that it started at -180. the delay always pushes
% the vectors counter clockwise. (This is simply mod(val, pi), MA20150727)
delay_hor  = delay_hor  + (pi/2)*(1-sign(delay_hor));
delay_vert = delay_vert + (pi/2)*(1-sign(delay_vert));

% angle-correct the delay maps according to parameter delaythr
delaythr = 0;
angind = delay_hor < delaythr;
delay_hor(angind) = abs(delay_hor(angind) - pi);
angind = delay_vert < delaythr;
delay_vert(angind) = abs(delay_vert(angind) - pi);

% use delay vector to calculate retinotopy.
map_hor   = 0.5*(angle(exp(1i*(ang0-delay_hor)))  - angle(exp(1i*(ang2-delay_hor))));
map_vert  = 0.5*(angle(exp(1i*(ang1-delay_vert))) - angle(exp(1i*(ang3-delay_vert))));

% radians to degrees
map_hor   = (map_hor  /(2*pi)) * visual_field(1);
map_vert  = (map_vert /(2*pi)) * visual_field(2);

%% Smooth the maps if needed and then get the visual field sign map
% GF = fspecial('gaussian',size(map_hor), 3);
% GF = GF/sum(GF(:));
% map_hor = ifft2( fft2(map_hor).*abs(fft2(GF)) );
% map_vert = ifft2( fft2(map_vert).*abs(fft2(GF)) );

[VFS, VFS_thr] = getVisSign(map_hor, map_vert);

%% Save the file for segmentation analysis
save([DIRS.SaveDir 'AnalyzedRet'],'CmplxMaps','EyeX', 'EyeY','visual_field')
% return;




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
for idx  = 1:length(Expts)
    iexp  = Expts(1,idx);
    for istim=1:p.nstim
        MyMin = prctile( VDAQ_ALL(iexp).repeatVar{istim}(:), 1);
        MyMax = prctile( VDAQ_ALL(iexp).repeatVar{istim}(:), 99);
        figure('color','w','WindowStyle', 'docked');
        for irep=1:VDAQ.nrepeats, 
           subplot(4,5,irep); imshow(fliplr(rot90(VDAQ_ALL(2).repeatVar{istim}(:,:,irep),3)),[MyMin MyMax]); axis image
        end
        %VarofVars = (var(VDAQ_ALL(2).repeatVar{1},[],3));
        %figure; imshow(fliplr(rot90(VarofVars,3)),[]); axis image
    end
end

%% superimpose on blood vessel patterns
%load([DIRS.SaveDir 'VisAreas_mouse15354.mat'])
%iexp=2;
dirPCO = [DIRS.camera '\', animal, '\' num2str(iexp) '\1\'];
pcoFiles  = dir( fullfile(dirPCO,'*.mat' ));
pcoFiles  = struct2cell(pcoFiles);
pcoFilesS = sort_nat(pcoFiles(1,:)); % sort files in natural order
pcoFilesS = pcoFilesS (1:numel(pcoFilesS)); % now ignore the last trials that are corrupted
filename = fullfile([dirPCO, pcoFilesS{1}]);
[nRows, nCols, timeStamps, rawData, startTime] = loadPCOFile(filename);
rawData = (fliplr(imrotate(rawData,-90)));
figure('color','w','WindowStyle', 'docked');
imshow(rawData(:,:,1),[]); axis image
hold on; imcontour(imresize(VFS_thr,[nCols nRows],'nearest'),[1,1],'r')
%hold on; imcontour(imresize(BW,[nCols nRows],'nearest'),[1,1],'r')

%% see the time traces of a specified roi
global loopcond; loopcond=1; 
ChooseROI_SeeTimeTrace(VDAQ.tensor(1),VFS_thr)





