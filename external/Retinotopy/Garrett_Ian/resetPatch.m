function [im splitflag Npatch] = resetPatch(im,centerPatch,imlab,q)

%First make a dilated version of the original patch, as defined by the
%Sereno 
idorig = find(imlab == q);
imorigpatch = zeros(size(im));
imorigpatch(idorig) = 1; %the original region
SE = strel('disk',1,0);
imdilpatch = imdilate(imorigpatch,SE);

%Now make the same patch, but limited to the pixels at the center of
%space
idpatch = find(imlab == q & centerPatch); %get pixels for this patch that are at the center of visual space
impatch = zeros(size(im));
impatch(idpatch) = 1; %an image of the patch we are looking at
SE = strel('disk',1,0);
impatch = imopen(impatch,SE);
idpatch = find(impatch);
imlabdum = bwlabel(impatch,4);
idlab  = unique(imlabdum);
Npatch = length(idlab)-1;

splitflag = 0;
if length(idlab) > 2 %Did limiting the patch to the center of v. space "split it"? i.e. should it be multiple areas?

    imdist = bwdist(impatch); %distance trx on "truncated patch
    id = find(~imdilpatch); %Make a boundary around the patch
    imdist(id) = -inf;
    imdist(idpatch) = 0; %force the local minima before watershed
    wshed = watershed(imdist,4);
    wshed = sign(phi(wshed-1));  %make it isolated patches

    SE = strel('disk',1,0);
    wshed = imerode(wshed,SE);  %erode slightly just so the fracture is a bit wider
    im(idorig) = 0;
    im = im+wshed;  %replace old patch with the "split" one

    splitflag = 1;
end
