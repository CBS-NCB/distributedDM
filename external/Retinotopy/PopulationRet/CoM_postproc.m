function cdat = CoM_postproc(cdat,fidx,fdat,pdat)
% take the individual session's centers of masses, translate/rotate them to
% align with a particular session (denoted by fidx), and export them in the
% same structure.
%
% inputs ------------------------------------------------------------------
%   cdat : #session x 1 structure, with raw session-by-session CoM data.
%          this should be either (1) the output of CoM_compute OR (2) the 
%          contents of the GUI output com.mat.
%   fidx : 1 x 1 integer, denoting which session of cdat to which all other
%          sessions will be aligned. if zero or empty, all sessions will be
%          aligned to the mean direction 
% outputs -----------------------------------------------------------------
%   cdat : #session x 1 structure, identical in composition to the input
%          but with its contents changed accordingly. V1's location will be
%          zero and all sessions will be rotated appropriately, s.t.
%          cdat(k).dir is equal for all k.
% -------------------------------------------------------------------------
% note: this function is extracted directly from the GUI.
% - MJM
disp('****Scroll through the output rotations manually and check for errors****');

% assign the angle to which each session will be aligned
if isempty(fidx) || fidx == 0
    rotdir = mean([cdat.dirV]);
else
    rotdir = cdat(fidx).dirV;
end
% draw the scatter-centers (after deleting the old ones)
for k = 1:numel(cdat)
    % find the translation
    transl = cdat(k).xy(cdat(k).V1idx,:);
    transl_idx = 0.*transl;
    transl_idx(1) = round((0.5-transl(1))/mean(unique(diff(fdat(k).xdom))));
    transl_idx(2) = round((0.5-transl(2))/mean(unique(diff(fdat(k).ydom))));
    % find the rotation
    rotang  = cdat(k).dir - rotdir;
%     if pdat(k).flipLR, rotang = rotang - pi/2; end;
    % translate and rotate the maps (note translate to [0.5 0.5] and keep
    % using xdom and ydom)
    maptmp = imtranslate(fdat(k).visap_us,transl_idx);
    cdat(k).map0 = imrotate(maptmp,rotang*180/pi,'nearest','crop');
        if ~isequal(size(cdat(k).map0),[178 159])
            cdat(k).map0 = imresize(cdat(k).map0,[178 159]);
        end
    vfstmp = imtranslate(fdat(k).VFS,transl_idx);
    cdat(k).vfs0 = imrotate(vfstmp,rotang*180/pi,'nearest','crop');
        if ~isequal(size(cdat(k).vfs0),[178 159])
            cdat(k).vfs0 = imresize(cdat(k).vfs0,[178 159]);
        end
    % translate and rotate the CoMs (note translate to [0 0] and use
    % arbitrary axis units)
    cdat(k).xy0 = bsxfun(@minus,cdat(k).xy,transl)*rotmat(rotang);
    cdat(k).dir0 = rotdir;
end