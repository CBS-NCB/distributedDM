function im = splitPatchesX(im,kmap_hor,kmap_vert,kmap_rad,pixpermm)

xsize = size(kmap_hor,2)/pixpermm;  %Size of ROI mm
ysize = size(kmap_hor,1)/pixpermm; 
xdum = linspace(0,xsize,size(kmap_hor,2)); 
ydum = linspace(0,ysize,size(kmap_hor,1)); 
[xdom, ydom] = meshgrid(xdum,ydum); %two-dimensional domain

kmap_rad = smoothPatchesX(kmap_rad,im); %smooth the larger patches

hh = fspecial('gaussian',size(kmap_hor),3);
kmap_horS  = ifft2(fft2(kmap_hor) .*abs(fft2(hh))); % smooth the maps of azimuth & altitude
kmap_vertS = ifft2(fft2(kmap_vert).*abs(fft2(hh)));

[dhdx, dhdy] = gradient(kmap_horS);
[dvdx, dvdy] = gradient(kmap_vertS);
Jac = (dhdx.*dvdy - dvdx.*dhdy)*pixpermm^2; %deg^2/mm^2  %magnific factor is determinant of Jacobian (MA: eq 4 Garrett 2014)

%% Make Interpolated data to construct the visual space representations%%%
dim = size(kmap_horS);
U = 3;
xdum = linspace(xdom(1,1),xdom(1,end),U*dim(2)); 
ydum = linspace(ydom(1,1),ydom(end,1),U*dim(1));
[xdomI, ydomI] = meshgrid(xdum,ydum); %upsample the domain
kmap_horI  = round(interp2(xdom,ydom,kmap_horS ,xdomI,ydomI,'spline'));
kmap_vertI = round(interp2(xdom,ydom,kmap_vertS,xdomI,ydomI,'spline'));
kmap_radI  = round(interp2(xdom,ydom,kmap_rad  ,xdomI,ydomI,'spline'));

[dhdx, dhdy] = gradient(kmap_horI);
[dvdx, dvdy] = gradient(kmap_vertI);
JacI = (dhdx.*dvdy - dvdx.*dhdy)*(pixpermm*U)^2; % MA: JacI is Jac interpolated to become U times larger

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
SE = strel('disk',1,0);
im = imopen(sign(im),SE);
imlab = bwlabel(im,4); 
imdom = unique(imlab);

SE = strel('disk',1,0);
im = imerode(sign(im),SE); % MA: erosion of each patch by a disk of r=1

sphdom = -90:90;  %create the domain for the sphere
[sphX, sphY] = meshgrid(sphdom,sphdom);

%% Find patches to 'split'
clear spCov
%R = 30; %Find local min within central 30 deg (increasing R may increase the number of local minima detected MA)

global extraParams;
R =extraParams.R;


%%%%First limit patches to a certain eccentricity, R
imlab = bwlabel(im,4); 
imdom = unique(imlab);


centerPatch = getCenterPatch(kmap_rad,im,R);
% load('Dum.mat')
% Dum2= imresize(Dum,39/157.8,'nearest');
% centerPatch = Dum2;
% centerPatchI  = imresize(Dum,size(xdomI));

for q = 1:length(imdom)-1 %loop through each patch ("visual area")    
    im = resetPatch(im,centerPatch,imlab,q);       
end

imI = round(interp2(xdom,ydom,im,xdomI,ydomI,'nearest')); %interpolate to be the same size as the maps

%%%%%Now proceed with splitting based on no. of minima
imlab = bwlabel(im,4); 
imdom = unique(imlab);

imlabI = bwlabel(imI,4); 

% centerPatch =  getCenterPatch(kmap_rad,im,R);
% centerPatchI = getCenterPatch(kmap_radI,imI,R);

for q = 1:length(imdom)-1 %loop through each patch ("visual area")
    
    idpatch = find(imlab == q & centerPatch); %if a patch doesn't have any centerpatch, idpatch is empty.        
    dumpatch = zeros(size(im));
    dumpatch(idpatch) = 1;
    
    idpatchI = find(imlabI == q & centerPatchI);   
    dumpatchI = zeros(size(imI));
    dumpatchI(idpatchI) = 1;

    
    %figure,imagesc(dumpatch)
    
    Nmin = 1;
    if ~isempty(find(idpatch))
        
        
        %%Determine if it has a overlapping representation of visual space%%%%%%

        [spCov, JacCoverage, ActualCoverage, MagFac] = overRep(kmap_horI,kmap_vertI,U,JacI,dumpatchI,sphdom,sphX,pixpermm);
        
        % MA: if JacCov/ActualCov ratio is larger than 1, the patch contains redundant representation
        CovOverlap = JacCoverage/ActualCoverage;

        
        if CovOverlap > .999 
            
            %figure, imagesc(dumpatch)
            
            hor_cent = median(kmap_hor(find(dumpatch)));
            vert_cent = median(kmap_vert(find(dumpatch)));
            
            kmap_rad_cent = sqrt((kmap_hor-hor_cent).^2 + (kmap_vert-vert_cent).^2);
            
            kmap_rad_dum = zeros(size(kmap_rad));
            kmap_rad_dum(idpatch) = kmap_rad_cent(idpatch);

            [Nmin, minpatch, centerPatch2, Rdiscrete] = getNlocalmin(idpatch,R,kmap_rad_dum);

            [im, splitflag, Nsplit] = resetPatch(im,centerPatch2,imlab,q);
        end

    end


    if Nmin > 1
        id = find(imlab == q);
        dumpatch = zeros(size(im)); dumpatch(id) = 1;

        figure,
        subplot(1,3,1), ploteccmapsplit(dumpatch.*kmap_rad_cent,[0 45],1,pixpermm);
        title('Smoothed eccentricity map'), colorbar off
        ylabel('mm'), xlabel('mm')

        id = find(Rdiscrete == median(Rdiscrete(:))); Rdiscrete(id) = 0;
        subplot(1,3,2), ploteccmapsplit(dumpatch.*Rdiscrete,[0 45],1,pixpermm);
        hold on, contour(xdom,ydom,minpatch,[.5 .5],'k')
        title(['Discretized map; ' num2str(Nmin) ' minima found']), colorbar off

        subplot(1,3,3), ploteccmapsplit(dumpatch.*kmap_rad_cent.*im,[0 45],1,pixpermm);
        title('Flood the patch with watershed')


    end


end

%% Compute level of over-representation

imlab = bwlabel(im,4); 
imdom = unique(imlab);

imI = round(interp2(xdom,ydom,im,xdomI,ydomI,'nearest')); %interpolate to be the same size as the maps
imI(find(isnan(imI))) = 0;
imlabI = bwlabel(imI,4);

% centerPatch = getCenterPatch(kmap_rad,im,R);
% centerPatchI = getCenterPatch(kmap_radI,imI,R);
clear spCov JacCoverage ActualCoverage MagFac

for q = 1:length(imdom)-1 %loop through each patch ("visual area")

    idpatch = find(imlab == q & centerPatch);
    dumpatch = zeros(size(im));
    dumpatch(idpatch) = 1;
    
    idpatchI = find(imlabI == q & centerPatchI);   
    dumpatchI = zeros(size(imI));
    dumpatchI(idpatchI) = 1;
    
    %figure,imagesc(dumpatchI)

    [spCov{q}, JacCoverage(q), ActualCoverage(q), MagFac(q)] = overRep(kmap_horI,kmap_vertI,U,JacI,dumpatchI,sphdom,sphX,pixpermm);

    CovOverlap = JacCoverage(q)/ActualCoverage(q);

    % if JacCoverage/(pi*R^2) < .01 %get rid of the areas with practically no screen coverage
    if CovOverlap > 1.05 || JacCoverage(q)/(pi*R^2) < .001

        id = find(imlab == q);
        im(id) = 0;
    end

end

figure,
scatter(JacCoverage,ActualCoverage)
hold on
plot([0 max(JacCoverage)], [0 max(JacCoverage)],'k')
xlabel('Jacobian integral (deg^2)')
ylabel('Actual Coverage (deg^2)')



%     if CovOverlap > 1.05 || JacCoverage(q)/(pi*R^2) < .01
% 
%         id = find(imlab == q);
%         im(id) = 0;
%     end


