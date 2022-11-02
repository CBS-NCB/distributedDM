%% load final area seg containing all the data
clear all; close all; clc
load('\\Labserver2\data\MOUSE\Segmentation\data\areasegres\final.mat')

%% plot for BTBD3
ReszF = 0.24715;
mmperpix = 158.7*ReszF;
clear Areas
N =numel(fdat);

%%
% figure
for i = 17
    subplot(4,6,i)
    % get axis in mm
    xdom = (0:size(fdat(i).VFS,2)-1)*mmperpix;
    ydom = (0:size(fdat(i).VFS,1)-1)*mmperpix;
    imagesc(xdom./1000,ydom./1000,fdat(i).VFS), axis image
    
    % do morphological processes on segmented areas
    Areas(i).BW = funct_morphAreaSeg(fdat(i).visap_us);
    hold on; imcontour(xdom./1000, ydom./1000, Areas(i).BW,'k');
    title(fdat(i).filename)
    BW = Areas(i).BW;
    %BW = imresize(BW,[672 640],'nearest');
    %save(['\\Labserver2\data\MOUSE\IMAGING\GCAMP\M160226_MA\ANALYZED\2\' 'VisAreas_mouse15355.mat'],'BW')
    %save(['\\Labserver2\data\MOUSE\Segmentation\PopBTBD3\Areas\' fdat(i).filename '.mat'], 'BW' )
end

%% Compute center of masses for BTBD3
cdat = [];
for i = 1:N
    % compute the centers of mass of the current areas
    [comxy,comsign,~] = func_getPatchCoM( Areas(i).BW, fdat(i).visap_s);
    % convert the CoMs from indices to mm (parameter-independent)
    comxy = [interp1(1:numel(fdat(i).xdom), fdat(i).xdom, comxy(:,1))...
             interp1(1:numel(fdat(i).ydom), fdat(i).ydom, comxy(:,2))];
    
    [V1id,~,V1map] = func_getV1loc(Areas(i).BW, comxy);
    % get V1 grad dir
    meangraddir_hor = mean(fdat(i).pgraddir_hor(V1map));
    tmp = struct('xy',comxy,'sign',comsign,'V1idx',V1id,'dir',meangraddir_hor,'filename',fdat(i).filename);
    cdat = cat(1, cdat, tmp);
end

%% align and plot BTBD3
refDat = 8;
markerref = {'o','x'};
figure
imagesc(fdat(refDat).xdom,fdat(refDat).ydom, (fdat(refDat).visap_s+1).*fdat(refDat).visap_us); 
colormap bone; axis image; hold on;

clear cdatAlgn
Algn = CoM_postproc(cdat,refDat);
plot(Algn(refDat).xy(Algn(refDat).sign,1),Algn(refDat).xy(Algn(refDat).sign,2),'ro','markersize',10,'markerfacecolor','r')
plot(Algn(refDat).xy(~Algn(refDat).sign,1),Algn(refDat).xy(~Algn(refDat).sign,2),'bs','markersize',10,'markerfacecolor','b')

iBTBD3 = 8:16;
for i=setdiff(iBTBD3,refDat)
    plot(Algn(i).xy0(Algn(i).sign,1)+Algn(refDat).xy(Algn(refDat).V1idx,1), ...
        Algn(i).xy0(Algn(i).sign,2)+Algn(refDat).xy(Algn(refDat).V1idx,2),'ro')
    plot(Algn(i).xy0(~Algn(i).sign,1)+Algn(refDat).xy(Algn(refDat).V1idx,1), ...
        Algn(i).xy0(~Algn(i).sign,2)+Algn(refDat).xy(Algn(refDat).V1idx,2),'bs')
end

%% align and plot non-BTBD3
refDat = 2;
markerref = {'o','x'};
figure
imagesc(fdat(refDat).xdom,fdat(refDat).ydom, (fdat(refDat).visap_s+1).*fdat(refDat).visap_us); 
colormap bone; axis image; hold on; % colormap parula; 

clear cdatAlgn
Algn = CoM_postproc(cdat,refDat);
plot(Algn(refDat).xy(Algn(refDat).sign,1),Algn(refDat).xy(Algn(refDat).sign,2),'ro','markersize',10,'markerfacecolor','r')
plot(Algn(refDat).xy(~Algn(refDat).sign,1),Algn(refDat).xy(~Algn(refDat).sign,2),'bs','markersize',10,'markerfacecolor','b')

iNonBTBD3 = [1:2 4:7 17:20];
for i=setdiff(iNonBTBD3,refDat)
    plot(Algn(i).xy0(Algn(i).sign,1)+Algn(refDat).xy(Algn(refDat).V1idx,1), ...
        Algn(i).xy0(Algn(i).sign,2)+Algn(refDat).xy(Algn(refDat).V1idx,2),'ro')
    plot(Algn(i).xy0(~Algn(i).sign,1)+Algn(refDat).xy(Algn(refDat).V1idx,1), ...
        Algn(i).xy0(~Algn(i).sign,2)+Algn(refDat).xy(Algn(refDat).V1idx,2),'bs')
end










%% draw the scatter-centers (after deleting the old ones)
% rotdir = mean([cdat.dir]); %align around the mean; %cdat(1).dir;
% for i = 1:numel(cdat)
%     % move V1 to zero, rotate to orient with session #fidx
%     cdat(i).xy = bsxfun(@minus,cdat(i).xy,cdat(i).xy(cdat(i).V1idx,:))*rotmat(cdat(i).dir - rotdir);
%     cdat(i).dir = rotdir;
% end

%%
% figure
% for i=1:numel(cdat)
%     plot(cdat(i).xy(:,1),cdat(i).xy(:,2),'ro','markerfacecolor','w')
%     hold on
% end


%%
% for k = 1:numel(cdat)
%     % move V1 to zero, rotate to orient with session #fidx
%     cdat(k).xy0 = bsxfun(@minus,cdat(k).xy,cdat(k).xy(cdat(k).V1idx,:)) * rotmat(cdat(k).dir - rotdir);
%     cdat(k).dir0 = rotdir;
% end












