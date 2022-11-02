function fileList  = MakePCOFileList(DIR,p)

global VDAQ

% analysisdir = fullfile(DIRS.camera,VDAQ.animal, 'VDAQtensor', num2str(iseries),  num2str(iexp));
rawdatadir  = DIR; %fullfile(DIRS.camera,VDAQ.animal,num2str(iseries),num2str(iexp));

if ~exist(rawdatadir,'dir')
    fprintf('Data directory %s does not (yet?) exist\n', rawdatadir);
    return
end

%% create fileList

DirContents    = dir(sprintf('%s/*.mat',rawdatadir));

% DirContents    = DirContents(3:end);
nfiles      = length(DirContents);

if nfiles == 0
    error('There are no files to load');
end

% -------- sort DirContents by date
datenums = zeros(nfiles,1);
for ifile = 1:nfiles
    datenums(ifile) = datenum(DirContents(ifile).date);
end
[sorteddates, sortorder] = sort(datenums);
DirContents = DirContents(sortorder);

% Load only what's available according to protocol
% nfiles = length(p.seqnums(:)); % AB commented it out March 16, 2011
DirContents = DirContents(1:nfiles);

irepeat = 1;
counter = []; n = 1;
StimDoneFlags = zeros(p.nstim,1);
fileList = cell(p.nstim,1); % a bit of initialization
for ifile = 1:nfiles
    if strfind(DirContents(ifile).name,VDAQ.animal) ~= 1
        error('File name does not start with VDAQ.animal name');
    end
    filename = DirContents(ifile).name;
    thenumbers = sscanf(filename, 'M%d_%d.mat'); % stimulus, 0, counter
    % added by ND 101207 to more generally parse the filename
    ttt = regexpi(filename, '(\d*)', 'match');
    thenumbers = [str2num(ttt{1});str2num(ttt{end})]; %modified by AP in order to allow to put a number in the name after date and initials
    
    if isfield( p, 'seqnums' )
       counter(n) = thenumbers(2);
       [istim, irepeat ] = find( p.seqnums == counter(n) ); n = n+1;
    else
       error('Cannot find field p.seqnum!!');
    end
    
    disp([istim, irepeat])% for debug ts 20110520
    fileList{istim,irepeat} = DirContents(ifile);
end

if length(unique(diff(sort(counter))))>1
    error('Some camera files are missing!!');
end

nAvailableRepeats = size(fileList,2); % can very well be <  expt.nrepeats

incompleterepeat = 0;
for istim = 1:p.nstim
    if isempty( fileList{istim, nAvailableRepeats} );
        incompleterepeat = 1;
    end
end
if incompleterepeat
    nAvailableRepeats = nAvailableRepeats-1;
end

fprintf('Created file list for %d stimuli, %d repeats\n',size(fileList,1),nAvailableRepeats);