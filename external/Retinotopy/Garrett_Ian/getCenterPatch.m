function centerPatch = getCenterPatch(kmap_rad,im,R)

id = find(kmap_rad<R);  %Find pixels near the center of visual space
centerPatch = zeros(size(im));
centerPatch(id) = 1;  %Make a mask for them
centerPatch = centerPatch.*im;  
SE = strel('disk',2,0);
centerPatch = imopen(centerPatch,SE); %clean it up
centerPatch = medfilt2(centerPatch,[3 3]); 
