% Visual Area Segmentation (a la Garret 2014 J Neurosci) 
% --creation    150507@ma
% --last update 160729@mjm
%
% load the tensors from trials of retinotopy experiments, and repackage
% them into (1) an organized structure of tensors and necessary params, and
% (2) complex maps combining horizontal and vertical drifting bars.
% 
% this new version replaces the old version, which relied on some nested
% scripts and the needlessly bulky TensorBuild from MatteoBox. instead,
% this uses loadWFtrials, the universal widefield file-loader. note that
% we did not eradicate the global variable issue.
%
% run the code in its entirety from the top. the code will prompt you for
% information as necessary.

clearvars; clear global

%% 01 // Declarations and setting required directories
global DIRS VDAQ
%SetBaseDir;

% select one of these manually
% MODE 1A -- standard retinotopy, not redoing anything already done
% params = struct('singletrials',false,'repFlag',false,'resFac',0.8,'LoCutFreq',0.05,'HiCutFreq',2,'spatfiltwid',3,'frame0list',1:2);
% MODE 1B -- standard retinotopy, overwriting any existing copies
params = struct('singletrials',false,'repFlag',false,'resFac',0.8,'LoCutFreq',0.05,'HiCutFreq',2,'spatfiltwid',3,'frame0list',1:2,'moco_flag', true);
    % AND THEN COMMENT OUT THE LINES FOLLOWED BY ************************
% MODE 2 -- single trial retinotopy
% params = struct('singletrials',false,'repFlag',false,'resFac',0.8,'LoCutFreq',0.05,'HiCutFreq',2,'spatfiltwid',3,'frame0list',1:2);
    % AND THEN COMMENT OUT THE LINES FOLLOWED BY ************************

params.guisave = true;
guiSaveDir = '\\Labserver5\data\MOUSE\Segmentation\guidata';
params.regsave = false;
regSaveDir = '\\Labserver5\data\MOUSE\Segmentation\animalseg';
% prompt user for session information (without hard coding it!)
prompt = {'Indicator:','Animal/date:','ID: (only if saving reg files)','Series number:','Experiments: (enter as "# spacebar #")','Server:'};
def    = {'GCAMP','M170802_RA','16296','1','1 2','5'};
answer = inputdlg(prompt,'Enter session info...',1,def);
[indicator,animal,idnum,iseries,Expts,Serv] = deal(answer{:});
    iseries = str2double(iseries); 
    Expts   = cellfun(@str2double,regexp(Expts,'[0-9]*','match'));
    Serv    = str2double(Serv);

func_DefaultDirs(indicator, animal, iseries, Serv);
CameraInfo = GetCameraInfo(animal, iseries);
% prompt user for which VDAQ to load, or to make a new one



%% 02 // load neural data
% [AbsMaps,AngleMaps,CmplxMaps] = deal(cell(size(Expts)));

for idx  = 1:length(Expts)
    
    canddir = fullfile(DIRS.camera,animal,'VDAQtensor',num2str(iseries),num2str(Expts(1,idx)));
    candfiles = dir(fullfile(canddir,'VDAQ*.mat')); candfiles = {candfiles.name}; candfiles = cat(2,{'[MAKE NEW VDAQ]'},candfiles);
    fileselect = listdlg('PromptString','Select VDAQ loading method:','SelectionMode','single','ListString',candfiles,'Name',canddir);
    if fileselect == 1
        loadfile = [];
        savefile = inputdlg('Define a savename for the new VDAQ file:','',1,{'VDAQ.mat'}); savefile = savefile{1};
    else
        loadfile = fullfile(canddir,candfiles{fileselect});
    end
    
    % get the experiment number
    iexp  = Expts(1,idx);
    fprintf('starting expt %d of %d...\n',idx,length(Expts));
    % load associated data (params)
    p = ProtocolLoad(animal,iseries,iexp); % ProtocolInspect(p);%M150428 %
    if isempty(p.blankstims), p.blankstims = p.nstim; end
    load([DIRS.VS '_' num2str(iexp)]);    % Load ScreenInfo
    % get data (neural images) filenames
    fileDir = fullfile(DIRS.camera,animal,num2str(iseries),num2str(iexp));
    triallist = dir(fullfile(fileDir,'*.mat')); triallist = {triallist.name};
    triallist = sort_nat(triallist);
    
    %% other experiments
    if strcmp('vmovie3sequentialGrating.x',p.xfile)
        p.pfiledurs = p.pars(1,:)/10;
    end
    
    %%
    if isempty(loadfile)
        %% 02 // load and populate VDAQ global variable
        % (if tensor doesn't already exist) initialize VDAQ
        fprintf('loading trials and building VDAQ...\n');
        VDAQ = struct('animal',animal,'iseries',iseries,'iexp',iexp,'Frame0List',...
            params.frame0list,'ResizeFactor',params.resFac,'durs',p.pfiledurs(1),...
            'nstim',p.nstim,'fileList',{triallist},'MeanIntensities',[],...
            'MmPerCameraPix',CameraInfo.MmPerCameraPix*CameraInfo.BotFocalLength/CameraInfo.TopFocalLength/params.resFac);

        % load and average tensors (this new block of code replaces Main_WideField
%         p.seqnums(:,6:end) = [];
        if ~params.singletrials % averaged across trials
            for k = 1:p.nstim
                fprintf('    stimulus %d of %d\n',k,p.nstim);
                W = loadWFtrials( fileDir, triallist(p.seqnums(k,:)),[], 'alltrials_flag',true, ...
                    'downsamp',params.resFac, 'PreFrm',params.frame0list, 'trialavg_flag',true, ...
                    'LoCutFreq',params.LoCutFreq, 'HiCutFreq',params.HiCutFreq, 'trimendframe',0,...
                    'moco_flag',false, 'filter_flag',false );
                VDAQ.tensor{1,k} = W.xavg;
                VDAQ.MeanIntensities(k,:) = W.bg';
            end
        else % single trials
            for k = 1:p.nstim
                for kk = 1:size(p.seqnums,2)
                    fprintf('    stimulus %d of %d\n',k,p.nstim);
                    W = loadWFtrials( fileDir, triallist(p.seqnums(k,kk)),[], 'alltrials_flag',true, ...
                        'downsamp',params.resFac, 'PreFrm',params.frame0list, 'trialavg_flag',false, ...
                        'LoCutFreq',params.LoCutFreq, 'HiCutFreq',params.HiCutFreq, 'trimendframe',0,...
                        'moco_flag',false, 'filter_flag',false );
                    VDAQ.tensor(kk,k) = W.x;
                    VDAQ.MeanIntensities(k,kk) = W.bg;
                end
            end
        end
            % make all tensor units the same size
            fprintf('VDAQ structure postprocessing...\n');
            minsz = cellfun(@size,VDAQ.tensor,repmat({3},size(VDAQ.tensor)));
            S = struct('type','()','subs',cellfun(@transpose,cellfun(@squeeze,mat2cell(cat(3,repmat({':'},size(minsz)),...
                repmat({':'},size(minsz)),repmat({1:min(minsz)},size(minsz))),ones(size(minsz,1),1),ones(size(minsz,2),1),3),'uni',0),'uni',0));
            VDAQ.tensor = cellfun(@subsref,VDAQ.tensor,num2cell(S),'uni',0);
        VDAQ.nsummedframes = [];
        VDAQ.nrepeats = size(p.seqnums,2);
        VDAQ.meanvalue = nanmean(VDAQ.MeanIntensities(:));
        VDAQ.durs    = p.pfiledurs(1);
        VDAQ.FrameRate = median(reshape(cellfun(@size,VDAQ.tensor,repmat({3},size(VDAQ.tensor)))/p.pfiledurs(1),[],1));
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
            meanresps = zeros(size(VDAQ.tensor));
            for istim = 1:p.nstim
                for iten = 1:size(VDAQ.tensor,1)
                    meanresps(iten,istim) = mean(VDAQ.tensor{iten,istim}(:));
                end
            end
            VDAQ.meanvalue = mean(meanresps(:));
        for istim = 1:p.nstim
            for iten = 1:size(VDAQ.tensor,1)
                VDAQ.tensor{iten,istim} = (VDAQ.tensor{iten,istim}/VDAQ.meanvalue)-1;
            end
        end
        % saving
        if ~params.singletrials % only save trial-averaged data
            fprintf('saving VDAQ...\n');
            mkdir(fullfile(DIRS.camera,animal,'VDAQtensor',num2str(iseries),num2str(iexp)));
            % save the VDAQ file
            save(fullfile(DIRS.camera,animal,'VDAQtensor',num2str(iseries),num2str(iexp),savefile),'VDAQ','p','-v7.3')
            % save a registration file for behavior (if flagged)
            if idx == 1 && params.regsave
                matches = dir(fullfile(regSaveDir,sprintf('R%s_*',idnum))); matches = {matches.name};
                matches = regexp(matches,'[A-Z0-9]*\_([0-9]*)','tokens');
                if ~isempty(matches)
                    ridx = max(cellfun(@str2double,cellfun(@cell2mat,vertcat(matches{:}),'uni',0)));
                else
                    ridx = 0;
                end
                sessregim = nanmean(VDAQ.tensor{1},3);
                save(fullfile(regSaveDir,sprintf('R%s_%d.mat',idnum,ridx+1)),'sessregim');
            end
        else
            fprintf('cannot save this version VDAQ -- too big...\n');
        end
    else
       fprintf('VDAQ save already detected! automatically loading that for you... overwrite manually if you want to reanalyze\n') 
       load(loadfile);
    end
    
    % PERFORM ALL MAJOR PROCESSING STEPS (frame trimming, frame0-blank
    % corrections, spatiotemporal filtering)
    [p,stimfreqs] = processRetVDAQ(p,params);
    % compute complex maps
    VDAQ_full = VDAQ;
    for k = 1:size(VDAQ_full.tensor,1)
        VDAQ.tensor = VDAQ_full.tensor(k,:);
        [AbsMaps{k,idx}, AngleMaps{k,idx}, CmplxMaps{k,idx}] = TensorFrequency([], stimfreqs ); % Get Azm and elv maps
    end
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

%% Save the file for segmentation analysis
fprintf('saving AnalyzedRet...\n');
save([DIRS.SaveDir 'AnalyzedRet'],'CmplxMaps','EyeX', 'EyeY','visual_field')
if params.guisave
    save(fullfile(guiSaveDir,sprintf('AnalyzedRet_%s-%s-%d',indicator,animal(1:end-3),iseries)),'CmplxMaps','EyeX', 'EyeY','visual_field')
end 
fprintf('... complete ...\n'); 

return;



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
trn = 1;
ang0 = angle(CmplxMaps{trn,1}{1});
ang2 = angle(CmplxMaps{trn,1}{2});
ang1 = angle(CmplxMaps{trn,2}{1});
ang3 = angle(CmplxMaps{trn,2}{2});

% phase elev/azim/sign maps (orig author: Ian Nauhaus)
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

% Smooth the maps if needed and then get the visual field sign map
% GF = fspecial('gaussian',size(map_hor), 3);
% GF = GF/sum(GF(:));
% map_hor = ifft2( fft2(map_hor).*abs(fft2(GF)) );
% map_vert = ifft2( fft2(map_vert).*abs(fft2(GF)) );

[VFS, VFS_thr] = getVisSign(map_hor, map_vert, params.resFac);

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

[VFS1, VFS_thr1] = getVisSign(map_hor_1, map_vert_1, params.resFac);

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
%{
% load([DIRS.SaveDir 'VisAreas_mouse15354.mat'])
dirPCO = [DIRS.camera '\', animal, '\' num2str(iseries) '\' num2str(iexp) '\'];
pcoFiles  = dir( fullfile(dirPCO,'*.mat' ));
pcoFiles  = struct2cell(pcoFiles);
pcoFilesS = sort_nat(pcoFiles(1,:)); % sort files in natural order
pcoFilesS = pcoFilesS (1:numel(pcoFilesS)); % now ignore the last trials that are corrupted
filename = fullfile([dirPCO, pcoFilesS{21}]);
[nRows, nCols, timeStamps, rawData, startTime] = loadPCOFile(filename);
rawData = (fliplr(imrotate(rawData,-90)));
sessregim = mean(rawData,3);
Dir = '\\Labserver5\data\MOUSE\Segmentation\animalseg\';
save([Dir 'R' idnum '_1'],'sessregim')
%}

figure; imshow(sessregim(:,:,1),[]); axis image
hold on; imcontour(imresize(VFS_thr,size(sessregim),'nearest'),[1,1],'r')
%hold on; imcontour(imresize(BW,size(sessregim),'nearest'),[1,1],'r')

%func_saveFig([],[],200,[],[idnum '_' animal '_' num2str(iseries) '_'],0,0,1)


%% check power spectra
func_tensorPowerSpectra(VDAQ.tensor{1},[],VDAQ.durs(1),1,'k',VFS_thr)


%%
% Anim = '16267';
% Dir = '\\Labserver5\data\MOUSE\Segmentation\animalseg\';
% load([Dir 'A' Anim '_1'])
% load([Dir 'R' Anim '_1'])
% figure; imshow(sessregim(:,:,1),[]); axis image
% hold on; imcontour(imresize(BW,size(sessregim),'nearest'),[1,1],'r')
% func_saveFig([],[],100,['Areas_' Anim],[],0,0,1)








