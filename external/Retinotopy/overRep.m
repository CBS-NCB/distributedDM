function [spCov, JacCoverage, ActualCoverage, MagFac] = overRep(kmap_hor,kmap_vert,U,Jac,patch,sphdom,sphX,pixpermm)

pixpermm = pixpermm*U;

N = length(sphdom);

posneg = sign(mean(Jac(find(patch))));
id = find(sign(Jac)~=posneg | Jac == 0);
Jac(id) = 0;
patch(id) = 0;
    
idpatch = find(patch);
JacCoverage = abs(sum(abs(Jac(idpatch))))/pixpermm^2; %deg^2

sphlocX = round(kmap_hor(idpatch));
sphlocX = sphlocX-sphdom(1)+1;
sphlocY = round(kmap_vert(idpatch));
sphlocY = sphlocY-sphdom(1)+1;
sphlocVec = N*(sphlocX-1) + sphlocY;

spCov = zeros(size(sphX)); %a matrix that represents the sphreen
spCov(sphlocVec) = 1;
spCov = imfill(spCov);
SE = strel('disk',1,0);
spCov = imclose(spCov,SE);
spCov = imfill(spCov);
%spCov = medfilt2(spCov,[3 3]);
ActualCoverage = sum(spCov(:)); %deg^2
MagFac = ActualCoverage/length(idpatch);