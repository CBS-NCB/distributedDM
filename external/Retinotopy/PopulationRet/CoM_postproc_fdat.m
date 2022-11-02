function fdat = CoM_postproc_fdat(fdat,cdat)
% take the individual session's centers of masses, translate/rotate them to
% align with a particular session (denoted by fidx), and export them in the
% same structure.
%
% inputs ------------------------------------------------------------------
%   cdat : #session x 1 structure, with raw session-by-session CoM data.
%          this should be either (1) the output of CoM_compute OR (2) the 
%          contents of the GUI output com.mat.
%   fdat :
% outputs -----------------------------------------------------------------
%   fdat : 
% -------------------------------------------------------------------------
% note: this function is extracted directly from the GUI.
% - MJM

% draw the scatter-centers (after deleting the old ones)
for k = 1:numel(cdat)
    % move V1 to zero, rotate to orient with session #fidx
    cdat(k).xy0 = bsxfun(@minus,cdat(k).xy,cdat(k).xy(cdat(k).V1idx,:))*rotmat(cdat(k).dir - rotdir);
    cdat(k).dir0 = rotdir;
end