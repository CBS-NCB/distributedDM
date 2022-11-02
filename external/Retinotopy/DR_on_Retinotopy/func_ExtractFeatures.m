function out = func_ExtractFeatures(CmplxMaps,visual_field,Threshold,pixpermm)

[map_hor,map_vert] = func_MakeAngleMaps(CmplxMaps,visual_field,Threshold);

[graddir_hor, graddir_vert, VFS] = func_getPatchGarrett(map_hor,map_vert,pixpermm);

[Jac, MagFact, Eccentricity, prefAxisMF, Distrtion] = func_GetMagFactors(map_hor,map_vert,pixpermm);

out = [map_hor(:) map_vert(:) graddir_hor(:) graddir_vert(:) VFS(:) ...
       Jac(:) MagFact(:) Eccentricity(:) prefAxisMF(:) Distrtion(:)];
