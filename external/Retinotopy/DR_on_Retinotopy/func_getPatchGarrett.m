function [graddir_hor, graddir_vert, VFS] = func_getPatchGarrett(kmap_hor,kmap_vert,pixpermm)

% Compute visual field sign map
scales = [3 2];
mmperpix = 1/pixpermm;

[dhdx, dhdy] = gradient(kmap_hor);
[dvdx, dvdy] = gradient(kmap_vert);

graddir_hor = atan2(dhdy,dhdx);
graddir_vert = atan2(dvdy,dvdx);

vdiff = exp(-1i*graddir_hor) .* exp(1i*graddir_vert); %Should be vert-hor, but the gradient in Matlab for y is opposite.
VFS = sin(angle(vdiff)); %Visual field sign map
id = find(isnan(VFS));
VFS(id) = 0;

hh = fspecial('gaussian',size(VFS),scales(1)); 
hh = hh/sum(hh(:));
VFS = ifft2(fft2(VFS).*abs(fft2(hh)));  %Important to smooth before thresholding below