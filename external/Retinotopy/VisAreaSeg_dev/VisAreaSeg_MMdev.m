% Visual Area Segmentation (a la Garret 2014 J Neurosci) 
% orig MA 20150507, 
% edit MM 20160722,
%
% this new version replaces the old version, which relied on some nested
% scripts and the needlessly bulky TensorBuild from MatteoBox. instead,
% this uses loadWFtrials, the universal widefield file-loader.
%
% no inputs or edits required, just run the whole script top-to-bottom, and
% you'll be prompted for info as necessary.
clearvars; clear global

%% Declarations and setting required directories
global DIRS VDAQ
SetBaseDir;

params = struct('singletrials',false,'repFlag',false,'resFac',0.8,'LoCutFreq',0.05,'HiCutFreq',2,'spatfiltwid',3,'frame0list',1:2);
% prompt user for session information (without hard coding it!)
prompt = {'Indicator:','Animal/date:','Series number:','Experiments: (enter as "# spacebar #")','Server:'};
def    = {'GCAMP','M160227_MA','1','1 2','2'};
answer = inputdlg(prompt,'Enter session info...',1,def);
[indicator,animal,iseries,Expts,Serv] = deal(answer{:});
    iseries = str2double(iseries); 
    Expts   = cellfun(@str2double,regexp(Expts,'[0-9]*','match'));
    Serv    = str2double(Serv);

func_DefaultDirs(indicator, animal, iseries, Serv);
CameraInfo = GetCameraInfo(animal, iseries);

%% here, Load the data and do basic prep using Main_WideField.m for all the iexps
% [AbsMaps,AngleMaps,CmplxMaps] = deal(cell(size(Expts)));

for idx  = 1:length(Expts)
    % get the experiment number
    iexp  = Expts(1,idx);
    fprintf('starting expt %d of %d...\n',idx,length(Expts));
    % load associated data (params)
    p = ProtocolLoad(animal,iseries,iexp); % ProtocolInspect(p);%M150428 %
    if isempty(p.blankstims), p.blankstims = p.nstim; end
    load([DIRS.ScreenInfoDir '_' num2str(iexp)]);    % Load ScreenInfo
    % get data filenames
    fileDir = fullfile(DIRS.camera,animal,num2str(iseries),num2str(iexp));
    triallist = dir(fullfile(fileDir,'*.mat')); triallist = {triallist.name};
    triallist = sort_nat(triallist);
    
%     if exist(fullfile(DIRS.camera,animal,'VDAQtensor',num2str(iseries),num2str(iexp),'VDAQ.mat'),'file') < eps
        % initialize VDAQ
        fprintf('loading trials and building VDAQ...\n');
        VDAQ = struct('animal',animal,'iseries',iseries,'iexp',iexp,'Frame0List',...
            params.frame0list,'ResizeFactor',params.resFac,'durs',p.pfiledurs(1),...
            'nstim',p.nstim,'fileList',{triallist},'MeanIntensities',[],...
            'MmPerCameraPix',CameraInfo.MmPerCameraPix*CameraInfo.BotFocalLength/CameraInfo.TopFocalLength/params.resFac);

        % load and average tensors (this new block of code replaces Main_WideField
        keyboard; p.seqnums(:,6:end) = [];
        for k = 1:p.nstim
            fprintf('    stimulus %d of %d\n',k,p.nstim);
            W = loadWFtrials( fileDir, triallist(p.seqnums(k,:)), 'alltrials_flag',false, ...
                'downsamp',params.resFac, 'PreFrm',params.frame0list, 'trialavg_flag',true, ...
                'LoCutFreq',params.LoCutFreq, 'HiCutFreq',params.HiCutFreq, 'trimendframe',0,...
                'moco_flag',false, 'filter_flag',false );
            VDAQ.tensor{1,k} = W.xavg;
            VDAQ.MeanIntensities(k,:) = W.bg';
        end

            % make all tensor units the same size
            fprintf('VDAQ structure postprocessing...\n');
            minsz = cellfun(@size,VDAQ.tensor,repmat({3},size(VDAQ.tensor)));
            S = struct('type','()','subs',mat2cell(cat(2,repmat({':'},numel(minsz),2),repmat({1:min(minsz)},numel(minsz),1)),ones(numel(minsz),1),3));
            VDAQ.tensor = cellfun(@subsref,VDAQ.tensor,num2cell(S)','uni',0);
        VDAQ.nsummedframes = [];
        VDAQ.nrepeats = size(p.seqnums,2);
        VDAQ.meanvalue = nanmean(VDAQ.MeanIntensities(:));
        VDAQ.durs    = p.pfiledurs(1);
        VDAQ.FrameRate = median(cellfun(@size,VDAQ.tensor,repmat({3},size(VDAQ.tensor)))/p.pfiledurs(1));
        % build VDAQ global to match with the rest of the code
        [nr,nc,nt]   = size(VDAQ.tensor{1});
        VDAQ.durs    = nt/VDAQ.FrameRate;
        VDAQ.nframes = nt;
        VDAQ.ny      = nr;
        VDAQ.nx      = nc;
        VDAQ.tt      = linspace(0, VDAQ.durs, nt);
        VDAQ.PCOdata_transposed = true;
        VDAQ.BuildDate = datestr(now);
        p.pars(1,:)  = ones(size(p.pars(1,:))).*(VDAQ.durs*10);
        % deltaF/F normalization
            % note this is for consistency, not because i think this is the
            % best method
            meanresps = zeros(p.nstim,1);
            for istim = 1:p.nstim
                meanresps(istim) = mean(VDAQ.tensor{istim}(:));
            end
            VDAQ.meanvalue = mean(meanresps);
        for istim = 1:p.nstim
            VDAQ.tensor{istim} = (VDAQ.tensor{istim}/VDAQ.meanvalue)-1;
        end
        % saving
        fprintf('saving VDAQ...\n');
        mkdir(fullfile(DIRS.camera,animal,'VDAQtensor',num2str(iseries),num2str(iexp)));
        save(fullfile(DIRS.camera,animal,'VDAQtensor',num2str(iseries),num2str(iexp),'VDAQ.mat'),'VDAQ','p','-v7.3')
%     else
%         fprintf('VDAQ save already detected! automatically loading that for you... overwrite manually if you want to reanalyze\n')
%         load(fullfile(DIRS.camera,animal,'VDAQtensor',num2str(iseries),num2str(iexp),'VDAQ.mat'));
%     end
    
        % trim invalid frames from tensor
        figure; meansig = cell(1,p.nstim);
        for istim = 1:p.nstim
            col = {'ko-'}; if istim == p.blankstims, col = {'ro-'}; end
            meansig{istim} = squeeze(mean(mean(VDAQ.tensor{istim})));
            plot(meansig{istim}, col{1}); hold on
        end
        title('Click to choose an appropriate Threshold (y) for trimming')
        [~, tshld2] = ginput(1); close;
        keepidx = cellfun(@gt,meansig,repmat({tshld2},size(meansig)),'uni',0);
        [~,tmpidx] = min(cellfun(@sum,keepidx));
        keepidx = repmat(keepidx(tmpidx),size(keepidx));
        S = struct('type','()','subs',mat2cell(cat(2,repmat({':'},numel(keepidx),2),keepidx'),ones(numel(keepidx),1),3));
        VDAQ.tensor = cellfun(@subsref,VDAQ.tensor,num2cell(S'),'uni',0);
        [nr,nc,nt]   = size(VDAQ.tensor{1});
        VDAQ.durs    = nt/VDAQ.FrameRate;
        VDAQ.nframes = nt;
        VDAQ.ny      = nr;
        VDAQ.nx      = nc;
        VDAQ.tt      = linspace(0, VDAQ.durs, nt);
        p.pars(1,:)  = ones(size(p.pars(1,:))).*(VDAQ.durs*10);
    % prepare p-file for filtering procedures
    if isempty(p.blankstims), p.blankstims = p.nstim; end
    p.pfilefreqs = p.pars(4,:)/100;
    stimfreqs = p.pfilefreqs;
    stimfreqs(p.blankstims) = [];
    % Temporal and Spatial Filtering
    func_TempFiltering(params.LoCutFreq, params.HiCutFreq);
    func_SpatFiltering(params.spatfiltwid);
    % Blank and Frame0 corrections
    blanklist = p.blankstims;
    func_CorrectBlank(blanklist, 1);
    func_CorrectFrm0();
    
    [AbsMaps{idx}, AngleMaps{idx}, CmplxMaps{idx}] = TensorFrequency([], stimfreqs ); % Get Azm and elv maps
%     VDAQ_ALL(idx)= VDAQ; % Keep both VDAQs in workspace in case you need it.
    
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
EyeY = -1; 
EyeX = 0;

%% Save the file for segmentation analysis
fprintf('saving AnalyzedRet...\n');
save([DIRS.SaveDir 'AnalyzedRet_debug'],'CmplxMaps','EyeX', 'EyeY','visual_field')
fprintf('... complete ...\n'); return;

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
%%
% load([DIRS.SaveDir 'VisAreas_mouse15354.mat'])
% dirPCO = [DIRS.camera '\', animal, '\' num2str(1) '\1\'];
% pcoFiles  = dir( fullfile(dirPCO,'*.mat' ));
% pcoFiles  = struct2cell(pcoFiles);
% pcoFilesS = sort_nat(pcoFiles(1,:)); % sort files in natural order
% pcoFilesS = pcoFilesS (1:numel(pcoFilesS)); % now ignore the last trials that are corrupted
% filename = fullfile([dirPCO, pcoFilesS{1}]);
% [nRows, nCols, timeStamps, rawData, startTime] = loadPCOFile(filename);
% rawData = (fliplr(imrotate(rawData,-90)));
% 
% figure; imshow(rawData(:,:,1),[]); axis image
% % hold on; imcontour(imresize(VFS_thr,[nCols nRows],'nearest'),[1,1],'r')
% hold on; imcontour(imresize(BW,[nCols nRows],'nearest'),[1,1],'r')
