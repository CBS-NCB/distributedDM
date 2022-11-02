% load('\\Labserver2\data\MOUSE\Segmentation\guidata\areasegres\res_GCAMP-M151125-1.mat'); % 15098
load('\\Labserver2\data\MOUSE\Segmentation\guidata\areasegres\res_GCAMP-M160215-1.mat'); % 15100
% window the vfs map
figure; imagesc(curr.VFS); colormap jet; colorbar; set(gca,'clim',[-1 1]);
hold on; contour(curr.visap_us,[1 1],'color','k','linewidth',2);
h=imrect;
ROI = round(h.getPosition);
% base plot
figure; ax1 = axes;
imagesc(curr.VFS(ROI(2):ROI(2)+ROI(4)-1,ROI(1):ROI(1)+ROI(3)-1)); axis image;
colormap jet; colorbar; set(gca,'clim',[-1 1],'xtick',[],'ytick',[]); freezeColors;
% plot areas on top
[~,bordr] = func_dilateareas(imresize(curr.visap_us(ROI(2):ROI(2)+ROI(4)-1,ROI(1):ROI(1)+ROI(3)-1),3));
bordr = imdilate(bordr,strel('diamond',1));
axes('color','none','position',get(ax1,'position'),'ydir','reverse','xtick',[],'ytick',[]);
hold on; [~,ch] = contour(bordr,[1 1],'color','k','linewidth',2); axis image;