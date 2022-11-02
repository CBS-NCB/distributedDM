function [Jac, MagFact, Eccentricity, prefAxisMF, Distrtion] = func_GetMagFactors(kmap_hor,kmap_vert,pixpermm) %#ok<*FNDEF>


hh = fspecial('gaussian',size(kmap_hor),3); 
hh = hh/sum(hh(:));
kmap_hor  = ifft2(fft2(kmap_hor) .*abs(fft2(hh)));
kmap_vert = ifft2(fft2(kmap_vert).*abs(fft2(hh)));

[dhdx, dhdy] = gradient(kmap_hor);
[dvdx, dvdy] = gradient(kmap_vert);
Jac = (dhdx.*dvdy - dvdx.*dhdy)*pixpermm^2; %deg^2/mm^2  %magnification factor is determinant of Jacobian (MA: eq 4 Garrett 2014)

vecH = dhdx + 1i*dhdy; 
vecV = dvdx + 1i*dvdy;
Res = abs(vecH).*exp(1i*(angle(vecH) + pi/2)*2) + abs(vecV).*exp(1i*(angle(vecV) + pi/2)*2);
Res = Res./( abs(vecH)+abs(vecV) );
Distrtion = abs(Res);
prefAxisMF = angle(Res)/2*180/pi;
id = find(prefAxisMF<0); 
prefAxisMF(id) = prefAxisMF(id)+180;

MagFact = 1./(sqrt(abs(Jac))*pixpermm);

az = (kmap_hor)*pi/180; %azimuth
alt = (kmap_vert)*pi/180; %altitude
Eccentricity = atan(sqrt(tan(az).^2 + (tan(alt).^2)./(cos(az).^2)))*180/pi;  %eccentricity
