%%
SE = strel('disk',3,0); 
for i = 1:5
    for j = 1:5
        Tmp = imerode(Areas{i,j},SE);
        AreasLoop(:,:,i*j) = Tmp;
    end
end
SumAreas = sum(AreasLoop,3);
SumAreas = SumAreas>9;

BW = bwlabel(SumAreas,4);

%% fignd the holes and fill
BW = bwlabel(Areas); Tmp = BW;
BW5 = imfill(BW==5,'holes');
BW(Tmp==5)=0;
BW = BW +BW5;

%%
imseg = BW>0;
Npad = 30;  %Arbitrary padding value.  May need more depending on image size and resolution
dim = size(imseg);
imsegpad = [zeros(dim(1),Npad) imseg zeros(dim(1),Npad)];
dim = size(imsegpad);
imsegpad = [zeros(Npad,dim(2)); imsegpad; zeros(Npad,dim(2))];

SE = strel('disk',10,0);
imbound = imclose(imsegpad,SE);

imbound = imfill(imbound); %often unnecessary, but sometimes there are small gaps need filling

% SE = strel('disk',5,0);
% imbound = imopen(imbound,SE);

SE = strel('disk',3,0);
imbound = imdilate(imbound,SE); %Dilate to account for original thresholding.
imbound = imfill(imbound);

%Remove the padding
imbound = imbound(Npad+1:end-Npad,Npad+1:end-Npad);
imbound(:,1) = 0; imbound(:,end) = 0; imbound(1,:) = 0;  imbound(end,:) = 0; 

%Only keep the "main" group of patches. Preveiously used opening (see above), but this is more robust:
bwlab = bwlabel(imbound,4);
labid = unique(bwlab);
for i = 1:length(labid)
   id = find(bwlab == labid(i));
   S(i) = length(id);
end
S(1) = 0; %To ignore the "surround patch"
[~, id] = max(S);
id = bwlab == labid(id);
imbound = 0*imbound;
imbound(id) = 1;
imseg = imseg.*imbound;

%This is important in case a patch reaches the edge... we want it to be smaller than imbound
imseg(:,1:2) = 0; imseg(:,end-1:end) = 0; imseg(1:2,:) = 0;  imseg(end-1:end,:) = 0; 


%% Morphological thinning to create borders that are one pixel wide

%Thinning
bordr = abs(imbound-imseg);
bordr = bwmorph(bordr,'thin',Inf);
bordr = bwmorph(bordr,'spur',4);

%Turn border map into patches
VisAreas = bwlabel(1-bordr,4);
VisAreas(find(VisAreas == 1)) = 0;
VisAreas = sign(VisAreas);
%%
figure
imshow(VisAreas,[])

%% superimpose VFS and Segmeted Areas 
figure
imagesc(VFS); axis image
hold on;
imcontour(VisAreas,'k')

%%
GF = fspecial('gaussian',size(VFS), 2);
GF = GF/sum(GF(:));

ReszF   = 0.24715;
AbsMap  = imresize(abs(CmplxMaps{1}{2}),ReszF);
AbsMapS = ifft2( fft2(AbsMap).*abs(fft2(GF)) );

%find the pixel distances from the imbound
[mx, my]= meshgrid(1:size(VFS,2), 1:size(VFS,1));
[ry,cx] = find(imbound == 1);
distMat = zeros(size(mx));
for rind = 1:size(mx,1)
    for cind = 1:size(mx,2)
        d = min( (ry-my(rind,cind)).^2 + (cx-mx(rind,cind)).^2);
        d = d(1);
        distMat(rind,cind) = d;
    end
end

figure; imagesc(imcomplement(mat2gray(distMat)).*VFSs.*AbsMapS); axis image;  colorbar; %caxis([-.01 .01]);
hold on; imcontour(VisAreas,'k');

% smooth things outside imbound
VFS_outbound = VFS;
GF = fspecial('gaussian',size(VFS), 5);
GF = GF/sum(GF(:));
VFS_outbound = ifft2( fft2(VFS_outbound).*abs(fft2(GF)) );
VFS(~imbound) = VFS_outbound(~imbound);

figure; imagesc(VFS); axis image; 
hold on; imcontour(VisAreas,'k');

mmperpix = 158.7*ReszF;
xdom = (0:size(VFS,2)-1)*mmperpix;
ydom = (0:size(VFS,1)-1)*mmperpix;

% axes in mm
imagesc(xdom./1000,ydom./1000,VFS), axis image
hold on; imcontour(xdom./1000,ydom./1000,VisAreas,'k');



%%
[nr, nc] = size(CmplxMapsRoi{1}{1});
VisAreas = imresize(VisAreas,[nr nc],'nearest');
VisAreas = bwlabel(VisAreas,4);
figure; imshow(VisAreas,[]);

%% plot area contours over delay maps
mmperpix = 1/157.8;
xdom = (0:size(map_hor,2)-1)*mmperpix;
ydom = (0:size(map_hor,1)-1)*mmperpix;


figure; subplot(2,2,1)
imagesc(xdom,ydom, map_hor.*MaskAmpId ,[-50 50]); axis image; colorbar; hold on,
contour(xdom,ydom, VisAreas,[.5 .5],'k'); title('delay map azm')

subplot(2,2,2); imagesc(xdom,ydom,map_vert.*MaskAmpId,[-50 50]); axis image; colorbar; hold on,
contour(xdom,ydom,VisAreas,[.5 .5],'k'); title('delay map alt')

subplot(2,2,3)
imagesc(xdom,ydom,delay_hor.*MaskAmpId,[0 pi]); axis image; colorbar; hold on,
contour(xdom,ydom,VisAreas,[.5 .5],'k'); title('delay map azm')

subplot(2,2,4); imagesc(xdom,ydom,delay_vert.*MaskAmpId,[0 pi]); axis image; colorbar; hold on,
contour(xdom,ydom,VisAreas,[.5 .5],'k'); title('delay map alt')

%% plot segmented areas on top of amp-phase maps
iMaps{1,1}  = CmplxMaps{1}{1}; 
iMaps{1,2}  = CmplxMaps{1}{2};
iMaps{2,2}  = CmplxMaps{2}{2};
iMaps{2,1}  = CmplxMaps{2}{1};
%close all; 
bpViewComplexMaps(iMaps); 
for i =[2 3 5 6]
    figure(3)
    subplot(2,3,i)
    hold on; [~,h] = imcontour(VisAreas,'r-'); set(h,'LineWidth',1)
end
%%
save([DIRS.SaveDir 'VisualAreas'],'VisAreas','SumAreas','Areas' );


















