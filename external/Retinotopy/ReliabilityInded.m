% Calculating SNR for each area given different stimulus used for retinotopy (Mohammad 20151201)
clear all; close all; clc; clear global

%% Declarations and setting required directories
global DIRS VDAQ
SetBaseDir;
% addpath(fullfile(DIRS.baseDir,'\LAB_CODE\MATLAB\MOUSE\WIDEFIELD\Main_WideField'));
indicator = 'GCAMP'; animal= 'M150923_MA'; iseries = 1; Expts = [1 2]; Serv = 1;
func_DefaultDirs(indicator, animal, iseries, Serv);

%% Load processed area seg (the population data of BTBD3)
Dir = DIRS.SaveDir; %'\\Labserver2\data\MOUSE\Segmentation\PopBTBD3\Areas\';
% Files = dir(fullfile(Dir,'*.mat'));
% Files = struct2cell(Files);
% FN = sort_nat(Files(1,:));
% filename = FN{9};
filename= 'VisAreas_mouse39.mat';

%% Do every thing here
for idx  = 1:length(Expts)
    % load the tensor and do initial processing
%     iexp  = Expts(1,idx);
%     RepFlag = 0;
%     resFac = .8;
%     Main_WideField;
%     VDAQ_ALL(idx)= VDAQ; % Keep both VDAQs in workspace in case you need it.
%     VDAQ=[];
    
    % load the areaseg
    load([Dir filename])
    BW = imresize(BW,[nr nc],'nearest'); % resize areas to tensor size
    BW = bwlabel(BW,4);
    
    % Calc Amp Spect
    for iSeg = 1:max(BW(:))
        Area = repmat(BW==iSeg,1,1,nt); % logical tensor of an area
        for iStim = 1:size(VDAQ_ALL(idx).tensor,2)
            %calc power spectrum for each area given different stim
            Tmp = VDAQ_ALL(idx).tensor{iStim}.*Area;
            Tmp(Tmp==0)=NaN;
            y{iSeg}(iStim,:)= squeeze(nanmean(nanmean(Tmp)));
            [Frq, Amp{idx,iSeg}(iStim,:)] = PowerTimeTrace(y{iSeg}(iStim,:), VDAQ_ALL(idx).durs);
            display(['[iSeg iStim]: ' num2str([iSeg iStim])])
            clear Tmp
        end
    end
    
    % calc SNR for each area given different stim
    for iSeg = 1:max(BW(:))
        for iStim = 1:size(VDAQ_ALL(idx).tensor,2)
            %Amp{idx,iSeg}(iStim,:) = Amp{idx,iSeg}(iStim,:) - Amp{idx,iSeg}(3,:); %subtract the blank spect
            idxFrq = find(Frq > 0.03);
            StimFrq = myfreqs(1);
            [~, id]= min(abs(Frq-StimFrq)); %id = find(Frq==StimFrq);
            stimAmp = Amp{idx,iSeg}(iStim, id);
            sideAmp = Amp{idx,iSeg}(iStim, [id-1 id+1]);
            SNR{idx}(iSeg,iStim) = stimAmp / mean(sideAmp);
            display(['[iSeg iStim]: ' num2str([iSeg iStim])])
        end
    end
    
    % Plot amps for each area
    figure
    col = {'k','b','r'};
    for iSeg = 1:max(BW(:))
        subplot(3,4,iSeg); hold on
        for iStim = 1:size(VDAQ_ALL(idx).tensor,2)
            h(iStim)= plot(Frq(idxFrq), Amp{idx,iSeg}(iStim, idxFrq), col{iStim}); hold on
            xlim([0 2]); xlabel('Freq (Hz)'); ylabel('Amplitude')
            set(gca,'tickdir','out'); set(gcf,'color','w'); box off
            plot([StimFrq StimFrq],ylim,'k:')
        end
        title(['Area ' num2str(iSeg) '; SNR: ' num2str(SNR{idx}(iSeg,1)) ', ' ...
            num2str(SNR{idx}(iSeg,2)) ', ' num2str(SNR{idx}(iSeg,3))])
        set(h,'linewidth',2)
        legend(h,'direct','reverse','blank'); legend boxoff
        
    end
    
    % Plot areas' labels
    imProperties = regionprops(BW,'centroid');
    reszF = 0.8; pixpermm = 158.7;
    mmperpix = 1/(pixpermm*reszF);
    ydom =(1:nr)*mmperpix;
    xdom =(1:nc)*mmperpix;
    subplot(3,4,max(BW(:))+1); cla; imagesc(xdom, ydom, BW);
    xlabel('mm'); ylabel('mm');
    set(gca,'xtick',1:4,'ytick',1:4,'xticklabel',1:4,'yticklabel',1:4)
    colormap gray; axis image;
    for Lbl =1:max(BW(:))
        Cntr = [round(imProperties(Lbl).Centroid(1)), round(imProperties(Lbl).Centroid(2))]*mmperpix;
        text(Cntr(1), Cntr(2), num2str(Lbl),'fontsize', 10, 'color','red')
    end
end

%% Power of the whole visual cortex
figure;
for idx  = 1:length(Expts)
    subplot(1,2,idx)
    Area = repmat(BW,1,1,nt); % logical tensor of an area
    for iStim = 1:size(VDAQ_ALL(idx).tensor,2)
        %calc power spectrum for each area given different stim
        Tmp = VDAQ_ALL(idx).tensor{iStim}.*Area;
        Tmp(Tmp==0)=NaN;
        yfull(iStim,:)= squeeze(nanmean(nanmean(Tmp)));
        [Frq, Ampfull(iStim,:)] = PowerTimeTrace(yfull(iStim,:), VDAQ_ALL(idx).durs);
        
        % calc SNR for each area given different stim
        idxFrq = find(Frq > 0.03);
        StimFrq = 0.5;
        [~, id]= min(abs(Frq-StimFrq)); %id = find(Frq==StimFrq);
        stimAmp = Ampfull(iStim, id);
        sideAmp = Ampfull(iStim, [id-1 id+1]);
        SNR_full(idx,iStim) = stimAmp / mean(sideAmp);
        
        clear Tmp
        h(iStim)= plot(Frq(idxFrq), Ampfull(iStim, idxFrq), col{iStim},'linewidth',2); hold on
        xlim([0 2]); xlabel('Freq (Hz)'); ylabel('Amplitude')
        set(gca,'tickdir','out'); set(gcf,'color','w'); box off
        plot([StimFrq StimFrq],ylim,'k:')        
        
        display(['iStim: ' num2str(iStim)])
    end
    title(['All Areas; SNR: ' num2str(SNR_full(idx,1)) ', ' ...
            num2str(SNR_full(idx,2)) ', ' num2str(SNR_full(idx,3))])
    legend(h,'direct','reverse','blank'); legend boxoff
end


