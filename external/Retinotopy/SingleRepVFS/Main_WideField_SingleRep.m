% Main Analyses for Widefield Imaging Data MA 20150507
% clear all; close all; clc; clear global
global VDAQ

%% Load Protocol and Params
p = ProtocolLoad(animal,iseries,iexp); % ProtocolInspect(p);%M150428 %

%% Remove corrupted 1st or last frames. See the plot before determining the tshld2
figure(101); clf; 
for istim = 1:p.nstim
    col = {'ko-'}; if istim == p.nstim, col = {'ro-'}; end
    plot(squeeze(mean(mean(VDAQ.tensor{istim}))), col{1}); hold on
    tshld2 = ((max(squeeze(mean(mean(VDAQ.tensor{istim})))) - min(squeeze(mean(mean(VDAQ.tensor{istim})))))/2)+ ...
        min(squeeze(mean(mean(VDAQ.tensor{istim}))));
    xl = xlim; plot([xl(1) xl(2)],[tshld2 tshld2],'b:')
    hold on
end


%% Trim the last corrupted frames (if <tshld2)
GoodFrm = [];
for istim = 1:p.nstim
    Idx = find( squeeze(mean(mean(VDAQ.tensor{istim}))) > tshld2 );
    GoodFrm(istim,:) = [Idx(1),Idx(end)];
end
[Frm_1, FrmEnd] = deal(max(GoodFrm(:,1)),min(GoodFrm(:,2)));
for istim = 1:p.nstim
    VDAQ.tensor{istim} = VDAQ.tensor{istim}(:,:,Frm_1:FrmEnd);    
end

%% Reset time parameters (if trimmed)
[nr,nc,nt]   = size(VDAQ.tensor{1});
VDAQ.durs    = nt/VDAQ.FrameRate;
VDAQ.nframes = nt;
VDAQ.ny      = nr;
VDAQ.nx      = nc;
VDAQ.tt      = linspace(0, VDAQ.durs, nt);
p.pars(1,:)  = ones(size(p.pars(1,:))).*(VDAQ.durs*10);
% tsec = linspace(0, nt/VDAQ.FrameRate, nt);
if isempty(p.blankstims), p.blankstims = p.nstim; end

%% Pars of frequency of stimulation
p.pfilefreqs = p.pars(4,:)/100;
myfreqs = p.pfilefreqs;
myfreqs(p.blankstims) = [];

% %% 
% %PlotFourietAmpSpect % plots amplitude spectrum of the ROI, need to choose ROI for each stim
% 
% %% Temporal and Spatial Filtering
% func_TempFiltering(0.05, 5); % define cutoff frequencies
% func_SpatFiltering(3); % define width of Gaussian filter
% 
% %% Blank and Frame0 corrections
% blanklist = p.blankstims;
% func_CorrectBlank(blanklist, 1)
% func_CorrectFrm0()



