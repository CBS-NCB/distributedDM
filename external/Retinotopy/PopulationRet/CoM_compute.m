function cdat = CoM_compute(varargin)
% compute the individual session's centers of masses. note: after this
% function, the analysis is not complete, i.e. the individual sessions
% cannot be compared/plotted against each other yet. see CoM_postproc
%
% inputs ------------------------------------------------------------------
% the user has two options:
% 1) nothing
% 2) fdat : 1 x #session structure, exported directly from the GUI in
%           final.mat and in the workspace upon "Write All". contents are
%           numerous, but necessary here are the visual area maps, stored 
%           in visap_us (binary image) and visap_s (signed +-1 according to
%           the area's sign.
% outputs -----------------------------------------------------------------
%   cdat : #session x 1 structure, with raw session-by-session CoM data.
% -------------------------------------------------------------------------
% note: this function is extracted directly from the GUI.
% - MJM

% if no data was supplied, find and load it
if isempty(varargin)
    [fnm,fpath] = uigetfile('*.mat','Select GUI output final.mat:','multiselect','off');
    load([fpath fnm]);
else
    fdat = varargin{1};
end
% session-by-session, compute the CoM information and accumulate it in cdat
cdat = [];
for k = 1:numel(fdat)
    % compute the centers of mass of the current areas
    [comxy,comsign,~] = func_getPatchCoM(fdat(k).visap_us,fdat(k).visap_s);
    % convert the CoMs from indices to mm (parameter-independent)
    comxy = [interp1(1:numel(fdat(k).xdom),fdat(k).xdom,comxy(:,1)) interp1(1:numel(fdat(k).ydom),fdat(k).ydom,comxy(:,2))];
    [V1id,~,V1map] = func_getV1loc(fdat(k).visap_us,comxy);
    % get V1 grad dir
    meangraddir_hor = mean(fdat(k).pgraddir_hor(V1map));
    meangraddir_vert = mean(fdat(k).pgraddir_vert(V1map));
    tmp = struct('xy',comxy,'sign',comsign,'V1idx',V1id,'dir',meangraddir_hor,'dirV',meangraddir_vert,'filename',fdat(k).filename);
    cdat = cat(1,cdat,tmp);
end