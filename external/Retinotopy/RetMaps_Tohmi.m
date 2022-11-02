% Retinotopic Maps (Kalatsky 20013) MA 20150507
clear all; close all; clc; clear global
%% Declarations
global DIRS VDAQ
indicator = 'YC';
DIRS = func_DefaultDirs(indicator);
animal  = 'M150401_MA'; iseries = 3; iexp    = 1; 

%% Run the initial analyses
Main_WideField;

%% Clip the tensor image to the ROI
Tensor = VDAQ.tensor;
Tensor = ClipROI(Tensor);

%% Get the time trace
swT = p.pars(2,1); 
Frms = swT+5 : swT+10;
figure
subplot(121)
for istim = 1:size (Tensor, 2)-1
    TimeCourse(istim,:) = mean(mean(Tensor{istim})) - mean(mean(Tensor{3})); hold on
    plot(1:nt,TimeCourse(istim,:)); hold on; 
    set(gca, 'xtick',[0 10 20 30 40 50],'xticklabel',[0 1 2 3 4 5],'xlim',[0 50])
    YL = ylim; plot([swT(1) swT(1)],[YL(1) YL(2)],'k:'); hold on
    xlabel('Time')
    pbaspect([2,1,1])
end
subplot(122)
plot(1:nt, mean(TimeCourse), 'r'); hold on
set(gca, 'xtick',[0 10 20 30 40 50],'xticklabel',[0 1 2 3 4 5],'xlim',[0 50])
YL1 = ylim; plot([swT swT],[YL1(1) YL1(2)],'k:'); hold on
xlabel('Time')
pbaspect([2,1,1])

%% Diff Map after blank correction
DiffMap1 = Tensor{1} - Tensor{2};
mymin = prctile(DiffMap1(:), 2.5);
mymax = prctile(DiffMap1(:), 97.5);
figure
imagesc(mean(DiffMap1(:,:,Frms),3)); 
axis image; colormap jet
title('Diff Map 1')

%% Save figure in the specified directory (Choose FigName Appropriately)
FigName = [{'TimeTrace'},{'PowSpect'}, {'DiffMap'}];
myDir   = ['\\Labserver\data\MOUSE\IMAGING\' indicator '\' animal '\VDAQtensor\' num2str(iseries) '\' num2str(iexp) '\Figs\']; 
if exist(myDir, 'dir'), display ('The directory already exists'); else display('The directory does not exist'); end
filename    = [animal '_' num2str(iseries) '_' num2str(iexp) '_' FigName{3}];
saveas(gcf, [myDir, filename], 'jpeg')

%%  Make Movie Using Spatially and Temporally Filtered VDAQ.tensor
close all; clc
Fig = figure('position', [100 200 1600 700]);
for istim   = 1: size(VDAQ.tensor,2);
    mymin   = prctile(VDAQ.tensor{istim}(:), 2.5);
    mymax   = prctile(VDAQ.tensor{istim}(:), 97.5);
    for frm = 1:nt  
        drawnow
        subplot(1, p.nstim, istim)
        imagesc(squeeze( VDAQ.tensor{1,istim}(:,:,frm) ), [mymin, mymax])
        title(['Time: ' num2str(frm*(1/floor(VDAQ.FrameRate))) ' Sec'],'fontsize',15)
        colormap('gray'); colorbar;
        axis image
        TheMovie(istim, frm) = getframe;
    end
end

%% Save the movies in the directory (labserver)
myDir       = ['\\Labserver\data\MOUSE\IMAGING\' indicator '\' animal '\VDAQtensor\' num2str(iseries) '\' num2str(iexp) '\']; 
for stimNum = 1: size(VDAQ.tensor,2)
    filename    = [animal '_' num2str(iseries) '_' num2str(iexp) '_' num2str(stimNum) '.avi'];
    movie2avi(TheMovie(stimNum,:), [myDir filename], 'compression','None', 'fps', floor(VDAQ.FrameRate));
    %winopen([myDir filename])
end





