

%%
close all; clc
Fig = figure('position', [100 200 1600 700]);
trialNum = 2;
for istim   = 1: size(VDAQ.tensor,2);
    SingleTrVDAQ = TensorBPTime(  0.05, 5, VDAQ.FrameRate, VDAQ.alldata{istim,trialNum}, [] );
    
    mymin   = prctile(SingleTrVDAQ(:), 2.5);
    mymax   = prctile(SingleTrVDAQ(:), 97.5);
    
    for frm = 1:nt  
        drawnow
        subplot(1, p.nstim, istim)
        imagesc(fliplr(imrotate(squeeze( SingleTrVDAQ(:,:,frm) ),-90)), [mymin, mymax])
        title(['Time: ' num2str(frm*(1/floor(VDAQ.FrameRate))) ' Sec'],'fontsize',15)
        colormap('gray'); colorbar;
        axis image
        TheMovie(istim, frm) = getframe;
    end
end

%% Save the movies in the directory (labserver)
for stimNum = 1: size(VDAQ.tensor,2)
    filename    = [animal '_' num2str(iseries) '_' num2str(iexp) '_' num2str(stimNum) '_' num2str(trialNum) '.avi'];
    movie2avi(TheMovie(stimNum,:), [DIRS.SaveDir filename], 'compression','None', 'fps', floor(VDAQ.FrameRate));
    %winopen([DIRS.SaveDir filename])
end
