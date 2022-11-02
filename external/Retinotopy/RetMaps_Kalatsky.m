% Retinotopic Maps (Kalatsky 20013) MA 20150507
clear all; close all; clc; clear global

global DIRS VDAQ
SetBaseDir;
% addpath(fullfile(DIRS.baseDir,'\LAB_CODE\MATLAB\MOUSE\WIDEFIELD\Main_WideField'));

%% Declarations
indicator = 'GCAMP'; animal  = 'M150512_MA_test'; iseries = 1; iexp = 2; 
func_DefaultDirs(indicator, animal, iseries);

%% Run the initial analyses
Main_WideField;

%% look at the maps of coherence
[AbsMaps, AngleMaps, FourierMaps] = TensorFrequency([], myfreqs );
bpViewComplexMaps( FourierMaps );

%% look at the retinotopic maps (azimuth or elevation)
Divs = {}; % division of the maps in complex space (performing subtraction in angular space)
Divs{1} = FourierMaps{1,1}./FourierMaps{2,1};
bpViewComplexMaps(Divs);

%%
Prods = {}; % Prodcuts of the maps in complex space (performing summation in angular space)
Prods{1} = FourierMaps{1,1}.*FourierMaps{2,1};
bpViewComplexMaps(Prods);

%% Save figure in the specified directory (Choose Names Appropriately)
Indicator = 'Emx1GCaMP3'; FigName = [{'PowerSpect'},{'PhaseMaps'}, {'Azm'}, {'Elv'}];
myDir = ['C:\Users\Mohammad-BL\Google Drive\MA\Mohammad\Figs\Retinotopy\', animal(1:8), Indicator, '\']; 
if exist(myDir, 'dir'), display ('The directory already exists'); else mkdir(myDir), end
saveas(gcf, [myDir, animal(1:8), num2str(iseries), '_', num2str(iexp), '_', FigName{1}], 'jpeg')

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
for stimNum = 1: size(VDAQ.tensor,2)
    filename    = [animal '_' num2str(iseries) '_' num2str(iexp) '_' num2str(stimNum) '.avi'];
    movie2avi(TheMovie(stimNum,:), [DIRS.SaveDir filename], 'compression','None', 'fps', floor(VDAQ.FrameRate));
    %winopen([DIRS.SaveDir filename])
end





