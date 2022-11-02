function ChooseROI_SeeTimT(im,patches,SwT)
%im is a cell containing 3D matrices (tensors)
if nargin<3, SwT=1; end

global loopcond;
if loopcond == false
    display('loopcondition must be true')
    clear loopcond
end

f= figure(101); set(f,'position',[300 300 800 400])
uicontrol('style','push','string','Exit','callback','loopcond = false;');

%%
if ~iscell(im)
    error('the input should be a cell containing the tensor')
elseif iscell(im)
    n= numel(im);
    subplot(1,2,1); imshow(std(im{1}, [],3),[]);
    if ~isempty(patches)
        hold on; [~,H] = imcontour(patches, 'r-'); set(H,'LineWidth',2)
    end
    h = imrect;
    % prevent rect to move beyond the figure axes
    addNewPositionCallback(h,@(p) []);%title(['roi: ' mat2str(p,3)])
    fcn = makeConstrainToRectFcn('imrect',get(gca,'XLim'),get(gca,'YLim'));
    setPositionConstraintFcn(h,fcn);
    
    %when true, plot the image and timetrace of the Roi
    Counter = 1;
    while loopcond
        pos = int32(round(getPosition(h)));
        for i = 1:n
            imRoi = im{i}(pos(2):pos(2)+pos(4)-1, pos(1):pos(1)+pos(3)-1,:);
            pause(0.5); subplot(n,2,i*2); drawnow
            plot(squeeze(mean(mean(imRoi))), 'ko-'); hold on; 
            subplot(n,2,i*2); YL = ylim;  
            plot([SwT SwT],[YL(1) YL(2)],'r--'); hold on 
            text(SwT,YL(1),'T0','color','red');
            xlim([1 size(im{i},3)]); pbaspect([3 1 1]);            
        end
        Counter = Counter+1; % every ten loop, clear plots
        if Counter==25
            subplot(1,2,1); title('press enter to clear time traces and continue', ...
                'fontsize',15); pause
            for i = 1:n, subplot(n,2,i*2); cla; Counter=1; end; 
        else subplot(1,2,1); title(['roi: ' mat2str(pos,3)],'fontsize',15);
        end
    end
end










