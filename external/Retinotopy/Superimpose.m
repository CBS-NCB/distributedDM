% Supperimpose Area Maps and Vein Patterns
global VDAQ

map_horz = Azm; 
map_vert = Elv;

%% Calculate visual field sign map VSgnMap
[hdx, hdy] = gradient(map_horz);
[vdx, vdy] = gradient(map_vert);

graddir_horz = atan2(hdy, hdx);
graddir_vert = atan2(vdy, vdx);

vdiff   = exp(1i*graddir_horz) .* exp(-1i*graddir_vert); % Should be vert-horz, but see the comment above.
VSgnMap = sin(angle(vdiff)); % Visual field sign map
id      = isnan(VSgnMap);
VSgnMap(id) = 0;

%% Spatial filtering (smoothing) the VSgnMap
Sigma = 10;
GF = fspecial('gaussian',size(VSgnMap), Sigma); % Gaussian Sigma = 10
GF = GF/sum(GF(:));
VSgnMap = ifft2( fft2(VSgnMap).*abs(fft2(GF)) );

%% Apply threshold to the smoothed VSgnMap to create discrete patches
AbsVFS  = abs(VSgnMap);
Thresh  = 1.5 * std(VSgnMap(:));
VS_Thr  = (sign(AbsVFS - Thresh/2) + 1)/2;  %threshold VSgnMap at (+-) 1.5 std

id      = find(VS_Thr);
imdum   = VS_Thr.*VSgnMap; 
imdum(id)= imdum(id)+1.1;

%% Superimpose
figure
Fused = imfuse(imdum, mean(VDAQ.tensor{1}(:,:,1:100),3),'diff');
hold on
imshow(Fused)

%% On a single trial tensor. See photobleaching first sections to load single trial tensor.
% Img = imresize(mean(Tens{1},3),0.8);
% figure(2); imshow(Img,[])
% figure(2); hold on; imcontour(VFS_thr,'r')
% pt = findobj(0,'type','patch');
% set(pt,'linewidth',2);





