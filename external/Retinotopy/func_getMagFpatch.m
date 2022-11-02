function MagFpatch = func_getMagFpatch(kmap_hor, kmap_vert, avrAmp, Prctile, pixpermm)

[JacIm, ~, ~] = getMagFactors(kmap_hor,kmap_vert,pixpermm);
GF      = fspecial('gaussian',size(JacIm), 5);
GF      = GF/sum(GF(:));
Jac     = ifft2( fft2(JacIm).*abs(fft2(GF)) );
MagFact = sqrt(1./abs(Jac));

%%
minAmp  = prctile(avrAmp(:),Prctile(1));
AmpMap  = zeros(size(avrAmp))*NaN;
Idx     = find(avrAmp>minAmp);
AmpMap(Idx) = avrAmp(Idx);
MaskAmpId   = ~isnan(AmpMap);

%%
mymin = prctile(MagFact(:), Prctile(2));
mymax = prctile(MagFact(:), Prctile(3));

Tmp         = MagFact; 
IdxL        = find(MagFact<mymin); 
Tmp(IdxL)   = NaN;
IdxU        = find(MagFact>mymax); 
Tmp(IdxU)   = NaN;

MagFpatch   = double(Tmp.*MaskAmpId);
MagFpatch   = MagFpatch>0;
MagFpatch   = imfill(MagFpatch,'holes'); 

%%
% [CoMxy, ~]  = getPatchCoM(Areas); % CoMxy, coordinates of Center of Mass
% V1id        = getV1id(Areas);
% 
% Vcnt(1) = map_hor(round(CoMxy(V1id,2)),round(CoMxy(V1id,1)));  %Get point in visual space at the center of V1
% Vcnt(2) = map_vert(round(CoMxy(V1id,2)),round(CoMxy(V1id,1)));
% az      = (map_hor - Vcnt(1))*pi/180; %azimuth in radian.  MA: eq5 Garrett 2014=> (H-<H>)
% alt     = (map_vert- Vcnt(2))*pi/180; %altitude in radian. MA: eq5 Garrett 2014=> (V-<V>)
% 
% map_rad = atan(  sqrt( tan(az).^2 + (tan(alt).^2)./(cos(az).^2)  )  )*180/pi; 
% 
% Ecc_MagF_Amp =avrAmp.*sqrt(1./abs(Jac)).*mat2gray(map_rad);