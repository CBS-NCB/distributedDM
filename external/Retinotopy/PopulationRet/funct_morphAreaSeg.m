function im = funct_morphAreaSeg(BW)

%% some morphological processing
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

%% Remove the padding
imbound = imbound(Npad+1:end-Npad,Npad+1:end-Npad);
imbound(:,1) = 0; imbound(:,end) = 0; imbound(1,:) = 0;  imbound(end,:) = 0; 

%% Only keep the "main" group of patches. Preveiously used opening (see above), but this is more robust:
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

% remove any potential pixels on the image boundary
imseg(:,1:2) = 0; 
imseg(:,end-1:end) = 0; 
imseg(1:2,:) = 0;  
imseg(end-1:end,:) = 0; 

% create one-pixel-wide borders between areas (morphological 'thinning')
bordr = abs(imbound-imseg);
bordr = bwmorph(bordr,'thin',Inf);
bordr = bwmorph(bordr,'spur',4);

% turn border map into patches
im = bwlabel(1-bordr,4);
im(find(im == 1)) = 0;
im = sign(im);




















