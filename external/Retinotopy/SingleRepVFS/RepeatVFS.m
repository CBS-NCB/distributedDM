% Visual Area Segmentation for each Repeat
clear all; close all; clc; clear global

%% Declarations and setting required directories
global DIRS VDAQ
SetBaseDir;
% addpath(fullfile(DIRS.baseDir,'\LAB_CODE\MATLAB\MOUSE\WIDEFIELD\Main_WideField'));
indicator = 'GCAMP'; animal= 'M151020_MA'; iseries = 2; Expts = [1 2]; Serv =2;
func_DefaultDirs(indicator, animal, iseries, Serv);
CameraInfo.HardwareBinning  = 1; % binning horz and vert
CameraInfo.MmPerCameraPix   = 0.0063 * CameraInfo.HardwareBinning;

%% screen info
for idx  = 1:length(Expts)
    iexp  = Expts(1,idx);
    p = ProtocolLoad(animal,iseries,iexp);
    load([DIRS.ScreenInfoDir '_' num2str(iexp)]); % Load ScreenInfo
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
EyeY = -15;
EyeX = -9;

%% File List for Azm and Alt
VDAQ.FrameRate  = 50;
VDAQ.nstim      = p.nstim;
VDAQ.animal     = animal;
VDAQ.iseries    = iseries;
VDAQ.Frame0List = [];
ResizeFactor    = 0.8;
VDAQ.ResizeFactor = ResizeFactor;

iexp = Expts(1);
pAzm = ProtocolLoad(animal,iseries,iexp); % ProtocolInspect(p);%M150428 %
azmDIR = fullfile(DIRS.camera,animal,num2str(iseries),num2str(iexp));
azmfileList     = MakePCOFileList(azmDIR,pAzm);
iexp = Expts(2);
pAlt = ProtocolLoad(animal,iseries,iexp); % ProtocolInspect(p);%M150428 %
altDIR = fullfile(DIRS.camera,animal,num2str(iseries),num2str(iexp));
altfileList     = MakePCOFileList(altDIR,pAlt);

%%
for irepAzm = 1:pAzm.nrepeats
    for irepAlt = 1:pAlt.nrepeats
        for idx = 1:length(Expts)
            iexp = Expts(idx);
            if iexp ==1, % VDAQ for azm
                VDAQ.iexp       = iexp;
                VDAQ.durs       = pAzm.pfiledurs(1);
                VDAQ.fileList   = azmfileList;
                % load and process Azm
                VDAQ.tensor = [];
                VDAQ.tensor = TensorLoadRepeat(irepAzm);
                Main_WideField_SingleRep;
                [~,~, C(irepAzm,irepAlt).CmplxMaps{iexp}] = TensorFrequency([], myfreqs );
            elseif iexp==2 % VDAQ for alt
                VDAQ.iexp = iexp;
                VDAQ.durs = pAlt.pfiledurs(1);
                VDAQ.fileList   = altfileList;
                % load and process Alt
                VDAQ.tensor = [];
                VDAQ.tensor = TensorLoadRepeat(irepAlt);
                Main_WideField_SingleRep;
                [~,~, C(irepAzm,irepAlt).CmplxMaps{iexp}] = TensorFrequency([], myfreqs );
            end           
        end
        display(['Azm Rep:' num2str(irepAzm) ', Alt Rep:' num2str(irepAlt)])
    end
end
%%
for irepAzm = 1:p.nrepeats
    for irepAlt = 1:p.nrepeats
        % Direct and Reverse maps for altitude and azimuth
        ang0 = angle(C(irepAzm,irepAlt).CmplxMaps{1}{1});
        ang2 = angle(C(irepAzm,irepAlt).CmplxMaps{1}{2});
        ang1 = angle(C(irepAzm,irepAlt).CmplxMaps{2}{1});
        ang3 = angle(C(irepAzm,irepAlt).CmplxMaps{2}{2});
        
        % Segmentation (Ian Nauhaus)
        %Find delay as the angle between the vectors
        delay_hor  = angle( exp(1i*ang0)+exp(1i*ang2) );
        delay_vert = angle( exp(1i*ang1)+exp(1i*ang3) );
        
        %The delay always pushes the vectors counter clockwise. (This is simply mod(val, pi), MA20150727)
        delay_hor  = delay_hor  + (pi/2)*(1-sign(delay_hor));
        delay_vert = delay_vert + (pi/2)*(1-sign(delay_vert));
        
        % additional corrections
        angind = delay_hor < pi/2;
        delay_hor(angind) = abs(delay_hor(angind) - pi);
        angind = delay_vert < pi/2;
        delay_vert(angind) = abs(delay_vert(angind) - pi);
        
        %Use delay vector to calculate retinotopy.
        map_hor  = 0.5*(angle(exp(1i*(ang0-delay_hor)))  - angle(exp(1i*(ang2-delay_hor))));
        map_vert = 0.5*(angle(exp(1i*(ang1-delay_vert))) - angle(exp(1i*(ang3-delay_vert))));
        
        map_hor  = (map_hor  /(2*pi)) * visual_field(1);
        map_vert = (map_vert /(2*pi)) * visual_field(2);
        
        % get the visual field sign map
        [VFS{irepAzm,irepAlt}, ~] = getVisSign(fliplr(rot90(map_hor,3)), fliplr(rot90(map_vert,3)));
        close
    end
end

%% plot only VFS maps
mmperpix=158.7*ResizeFactor;
xdom = (0:size(VFS{1},2)-1)*mmperpix;
ydom = (0:size(VFS{1},1)-1)*mmperpix;
for irep= 1: p.nrepeats*p.nrepeats
    sp= mod(irep-1,p.nrepeats)+1;
    if sp==1, figure; end
    
    drawnow
    subplot(4,5,sp)
    imagesc(xdom, ydom, VFS{irep}); axis image
    display(['VFS for Repeat' num2str(irep)])
    
    allVFS(:,:,irep)=VFS{irep};
end

%% vfs z-score and std
figure
subplot(221); imagesc(mean(allVFS,3)), axis image; colorbar; title('mean VFS')
subplot(222); imagesc(std(allVFS,[],3)), axis image; colorbar; title('sd VFS')
zScore = (allVFS - mean(allVFS(:)))/std(allVFS(:));
subplot(223);imagesc(mean(allVFS,3)./std(allVFS,[],3)), axis image; colorbar; title('mu/sd')
subplot(224);imagesc(mean(zScore,3)), axis image; colorbar; title('z VFS')

%%
figure, subplot(121)
imagesc(abs( mean(allVFS,3)+1i*std(allVFS,[],3) )); axis image
title('abs(mu+1i.sd)')
subplot(122)
imagesc(angle( mean(allVFS,3)+1i*std(allVFS,[],3) )); axis image
title('angle(mu+1i.sd)')

%% use VFS calcualted from trial averaged tensor (need to load the analyzedret.mat)
% VFSmu = VFS1;
% figure, subplot(121)
% imagesc(abs( VFSmu+1i*std(allVFS,[],3) )); axis image; colorbar
% title('abs(mu+1i.sd)')
% subplot(122)
% imagesc(angle( VFSmu+1i*(std(allVFS,[],3)) )); axis image; colorbar
% title('angle(mu+1i.sd)')
% 
% NewAreaSeg= abs( VFSmu+1i*std(allVFS,[],3) )>0.3;
% figure; imshow(NewAreaSeg);

%%
% Movie of the trial
FrameRate =5;
figure; clf
mymin   = prctile(zScore(:), 2.5);
mymax   = prctile(zScore(:), 97.5);
for frm = 1:size(zScore,3)
    drawnow
    imagesc(zScore(:,:,frm), [mymin, mymax])
    axis image
    TheMovie(1,frm) = getframe;
end
close;
% Save the movies in the directory (labserver)
filename    = [animal '.avi'];
if exist([DIRS.SaveDir 'VFS_Rep\'], 'dir'), display ('Save Directory Already Exists');
else display ('No Save Directory. mkdir ...'); mkdir([DIRS.SaveDir 'VFS_Rep\']), end
movie2avi(TheMovie(1,:), [DIRS.SaveDir 'VFS_Rep\' filename], 'compression','None', 'fps', FrameRate);

%%
tic
save([DIRS.SaveDir 'VFS_Rep\VFS_Rep.mat'], 'VFS', 'C','-v7.3')
toc







