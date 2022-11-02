function scatplotSM(data,mdir,paramstruct) 
% SCATPLOTSM, SCATterPLOT, 
%   Steve Marron's matlab function
%     Makes Scatterplot Matrices (or a single Scatterplot),
%     based on a given set of direction vectors
%     This puts 1-d projections on the diagonals
%         (ordering same as direction vectors)
%     And corresponding pairwise projection scatterplots
%         off of the diagonals
%     Note: for non-orthogonal directon vectors, this shows projections
%               onto the two dimensional plane determined by the direction
%               vectors.  The first direction determines the horizontal axis,
%               but the vertical axis is orthogonal, which can generate
%               some non-intuitive visualizations (especially exaggerated
%               when the directions are close to each other).  This is the
%               default because it gives the best views for matrices of
%               scatterplots.  To overide this, and force a naive scatterplot
%               of the two 1-d projections, use    iforcenaivesp = 1 
%
% Inputs:
%   data    - d x n matrix of data, each column vector is 
%                    a "d-dim digitized curve"
%
%   mdir    - a set of d dimensional "direction vectors"
%                 which define axes of scatterplots
%                 (need not be orthogonal)
%                 number of columns (+ abs(npcadiradd)) determines 
%                 size of output scatterplot matrix,
%                 with 1d projections on the diagonal,
%                 and 2d projections off the diagonal
%                   Number of columns + abs(npcadiradd) must be >= 2
%                   (otherwise use projplot1SM.m or projplot2SM.m)
%                 mdir can be empty (but then should set npcadiradd)
%
%   paramstruct - a Matlab structure of input parameters
%                    Use: "help struct" and "help datatypes" to
%                         learn about these.
%                    Create one, using commands of the form:
%
%       paramstruct = struct('field1',values1, ...
%                            'field2',values2, ...
%                            'field3',values3) ;
%
%                          where any of the following can be used,
%                          these are optional, misspecified values
%                          revert to defaults
%
%                    Version for easy copying and modification:
%     paramstruct = struct('',, ...
%                          '',, ...
%                          '',) ;
%
%    fields            values
%
%    npcadiradd       Number of Principal Component Directions to Add
%                     to mdir.
%                     0  (default)  Don't add any, just directions in mdir
%                     integer > 0   Add this number of PC directions to mdir
%                         Note this assumes that no elements of mdir are PCA
%                         directions.  Otherwise will get an error.
%                     integer < 0   Add this number of PC directions,
%                         for the projections on the subspace orthogonal
%                         to the subspace generated by mdir
%
%    irecenter        0  leave data as is, no recentering
%                     1  (default)  do mean recentering
%                             Generally recommended, to 
%                               make axes appear sensible
%
%    iorient          orientation of added pc directions:
%                     0  use orientation (+ or -) that comes from Matlab
%                             eigen decomposition
%                     1  re-orient each eigenvector (i.e. multiply
%                             by +-1), to make max projection go in positive 
%                             direction.
%                     2  re-orient each eigenvector to make largest 2nd moment
%                             (on each side of the origin) go in positive direction
%                     3  (default)  re-orient each eigenvector to point mostly 
%                             towards the "45 degree line", 
%                             i. e. so that the sum of the entries is positive
%                             Note: this tends to keep data scores
%                             "mostly in the same quadrant" making it
%                             useful for non mean centered data views
%                     Note:  this only effects added pc directions 
%                             (not input mdir directions)
%
%    iforcenaivesp    Indicator to FORCE NAIVE ScatterPlot (instead of default
%                     of projection onto 2-d plane determined by two direction
%                     vectors)
%                     0  use default of projection on 2-d plane
%                     1  force naive scatterplot, with
%                          projections on 1st direction on horizontal axis
%                          and projections on 2nd direction on vertical axis
%                              (thus NOT a 2-d plane projection)
%                     This only has effect when input direction vectors
%                     are not orthogonal
%
%    icolor           0  fully black and white version (everywhere)
%                     string (any of 'r', 'g', 'b', etc.) that single color
%                     1  (default)  color version (Matlab 7 color default)
%                     2  time series version (ordered spectrum of colors)
%                     nx3 color matrix:  a color label for each data point
%                             (to be used everywhere, except SiZer & QQ
%                              useful for comparing classes)
%
%    markerstr        Can be either a single string with symbol to use for marker,
%                         e.g. 'o' (default), '.', '+', 'x'
%                         (see "help plot" for a full list)
%                     Or a character array (n x 1), of these symbols,
%                         One for each data vector, created using:  strvcat
%
%    isubpopkde       0  (default) construct kde using only the full data set
%                     1  partition data into subpopulations, using the color
%                            indicators in icolor (defaults to 0, unless icolor
%                            is an nx3 color matrix), as markers of subsets.
%                            The corresponding mixture colors are then used in
%                            the subdensity plot, and overlaid with the full 
%                            density shown in black
%                     2  Show only the component densities (in corresponding 
%                            colors), without showing the full population
%                            density
%
%    idataconn        indices of data points to connect with line segments
%                     []  (default) for not connecting any data points
%                     otherwise : x 2 matrix of pairs of indices of data points
%                     (thus all intergers from 1,...,n).
%                     For time series data, this can give a clearer view of the 
%                     data ordering by using [[1, 2];[2, 3];[3, 4];...].
%                     For bias adjustment, with matched pairs, each row should
%                     have the ind0ces of the matches.
%
%    idataconncolor   can be any of 'r', 'g', 'b', etc., to use that color for all
%                     default is 'k'
%                     or can be 2 for easy rainbow coloring, 
%                         intended for time series of curves
%                         (Caution: this will use the first part of icolor,
%                          so might make most sense to use with icolor = 2, 
%                          to avoid strange results)
%                     or can be color matrix, where the number of rows  
%                     is the same as the number of rows of idataconn
%                         (has no effect for idataconn = [])
%
%    idataconntype    can be any of '-', '--', '-.', ':'
%                     default is '-'
%                     or can be character array (created using strvcat) of these, 
%                     where the number of rows is the same as 
%                     the number of rows of idataconn
%                         (has no effect for idataconn = [])
%
%    ibigdot          0  (default)  use Matlab default for dot sizes
%                     1  force large dot size in prints (useful since some
%                              postscript graphics leave dots too small)
%                              (Caution: shows up as small in Matlab view)
%                              Only has effect when markerstr = '.' 
%
%    idatovlay        0  Do not overlay data on kde plots (on diagonal)
%                     1  (default) overlay data using heights based on data ordering
%                              Note:  To see "c.d.f. style" increasing line, 
%                                     should also sort the data
%                     2  overlay data using random heights
%                     another integer > 0,  overlay data, using random heights,
%                                           with this numbers as the seed (so can 
%                                           better match data points across plots),
%                                           (should be an integer with <= 8 digits)
%
%    ndatovlay     number of data points overlayed (only has effect for idatovlay > 0)
%                       1  -  overlay up to 1000 points 
%                                           (random choice, when more)
%                       2  -  (default) overlay full data set
%                       n > 2   -  overlay n random points
%
%    datovlaymax      maximum (on [0,1] scale, with 0 at bottom, 1 at top of plot)
%                     of vertical range for overlaid data.  Default = 0.6
%
%    datovlaymin      minimum (on [0,1] scale, with 0 at bottom, 1 at top of plot)
%                     of vertical range for overlaid data.  Default = 0.5
%
%    legendcellstr    cell array of strings for legend (nl of them),
%                     useful for (colored) classes, create this using
%                     cellstr, or {{string1 string2 ...}}
%                         Note:  These strange double brackets seems to be needed
%                                for correct pass to subroutine
%                                It may change in later versions of Matlab
%                     CAUTION:  If are updating this field, using a command like:
%                         paramstruct = setfield(paramstruct,'legendcellstr',...
%                     Then should only use single braces in the definition of
%                     legendecellstr, i. e. {string1 string2 ...}
%                     Also a way to add a "title" to the first plot
%                             for this, use input of form:  {{string}}
%                     Also can indicate symbols, by just adding (at least 
%                             for +,x.o) into the text
%                     Note:  this only appears on the first plot
%
%    mlegendcolor     nl x 3 color matrix, corresponding to cell legends above
%                     (not needed when legendcellstr not specified)
%                     (defaults to black when not specified)
%
%    maxlim        Matrix of axis limits
%                        Use [] for default of all automatically chosen, by axisSM
%                        Use 1 for symmetrically chosen, by axisSM
%                            (often preferred for centered plots, as in PCA)
%                        Otherwise, must be (size(mdir,2) + npcadir) x 2 
%                            matrix of axis limits, 
%                            with each row corresponding to a direction
%                        Note:  Use is generally not recommended,
%                        because defaults give "good visual impression
%                        of decomposition".  It is mostly intended to allow
%                        the highlighting of "visually different scales" in data.
%                        But this comes at the cost of reduced detail being visible in the plots
%
%    iplotaxes        0 do not plot axes
%                     1 (default) plot axes as determined by direction vectors, 
%                           using thin line type
%
%    iplotdirvec      0 (default) do not plot direction vectors
%                     1 plot direction vectors, using thick line type
%
%    ibelowdiag       0 leave off scatterplots below diagonal
%                     1 (default) show scatterplots both above and below diagonal
%
%    titlecellstr     cell array for making subplot titles
%                     default is an empty cell array, {} for no titles
%                     To add titles, this should be a cell array of vertical
%                     cellstrs, where each cell location corresponds to a subplot
%                     Can create this using a command like:
%                               {{strvcat('Title Top Left','Title Bottom Left') ...
%                                 strvcat('Title Top Right','Title Bottom Right')}}
%                          Careful: note transpose structure
%                     Can create a single row of titles accross the top with:
%                               {{'title1' 'title2'}}
%                     To skip titles on some plots, put a space string ' '
%                     in those locations
%
%    titlefontsize    font size for title
%                           (only has effect when plot is made here,
%                            and when the titlecellstr is nonempty)
%                     default is empty, [], for Matlab default
%
%    labelcellstr    Vertical cell array of strings for axis labels
%                        create this using cellstr, 
%                        or {{string1; string2; ...}}
%                            Note:  These strange double brackets seems to be 
%                                needed for correct pass to subroutine
%                                It may change in later versions of Matlab
%                    default is an empty cell array, {}, which then gives
%                        "Direction i", for the ith Direction Vector
%                    For no labels, use {''; ''; ''}
%                    Number of rows should be size(mdir,2) + abs(npcadiradd)
%                        Empty entries default to "Direction i"
%                        Left over entries are ignored
%                    Note labels correspond to columns of mdir
%                        (not the usual x and x labels)
%
%    labelfontsize    font size for axis labels
%                                    (only has effect when plot is made here,
%                                     and when a label str is nonempty)
%                           default is empty [], for Matlab default
%
%    savestr          string controlling saving of output,
%                         either a full path, or a file prefix to
%                         save in matlab's current directory
%                         Will add .ps, and save as either
%                             color postscript (icolor ~= 0)
%                         or
%                             black&white postscript (when icolor = 0)
%                         unspecified:  results only appear on screen
%
%    iscreenwrite     0  (default)  no screen writes
%                     1  write to screen to show progress
%
%         These parameters create data for use by Marc Niethammer's mexplorer
%         They require savestr, and either icolor or markerstr to be manually set
%         Also, all but cellsubtypes must be non-empty
%         Then they create a figure file, which is used by mexplorer
%
%    celltypes 
%    cellsubtypes
%    slidenames
%    slideids 
%
%
%
% Outputs:
%     Graphics in current Figure
%     When savestr exists,
%        Postscript files saved in 'savestr'.ps
%                 (color postscript for icolor ~= 0)
%                 (B & W postscript for icolor = 0)
%
% Assumes path can find personal functions:
%    bwsjpiSM.m
%    kdeSM.m
%    lbinrSM.m
%    vec2matSM.m
%    pcaSM.m
%    projplot1SM.m
%    projplot2SM.m
%    bwrfphSM.m
%    bwosSM.m
%    rootfSM
%    bwrotSM.m
%    bwsnrSM.m
%    iqrSM.m
%    cquantSM.m
%    axisSM.m


%    Copyright (c) J. S. Marron 2004-2013



%  First set all parameters to defaults
%
npcadiradd = 0 ;
irecenter = 1 ;
iorient = 3 ;
iforcenaivesp = 0 ;
icolor = 1 ;
markerstr = 'o' ;
isubpopkde = 0 ;
idataconn = [] ;
idataconncolor = 'k' ;
idataconntype = '-' ;
ibigdot = 0 ;
idatovlay = 1 ;
ndatovlay = 2 ;
datovlaymax = 0.6 ;
datovlaymin = 0.5 ;
legendcellstr = {} ;
mlegendcolor = [] ;
maxlim = [] ;
iplotaxes = 1 ;
iplotdirvec = 0 ;
ibelowdiag = 1 ;
titlecellstr = {} ;
titlefontsize = [] ;
labelcellstr = {} ;
labelfontsize = [] ;
savestr = [] ;
iscreenwrite = 0 ;

celltypes = [];
cellsubtypes = [];
slidenames = [];
slideids = [];

%  Now update parameters as specified,
%  by parameter structure (if it is used)
%
if nargin > 2 ;   %  then paramstruct is an argument

  if isfield(paramstruct,'npcadiradd') ;    %  then change to input value
    npcadiradd = getfield(paramstruct,'npcadiradd') ; 
  end ;

  if isfield(paramstruct,'irecenter') ;    %  then change to input value
    irecenter = getfield(paramstruct,'irecenter') ; 
  end ;

  if isfield(paramstruct,'iorient') ;    %  then change to input value
    iorient = getfield(paramstruct,'iorient') ; 
  end ;

  if isfield(paramstruct,'iforcenaivesp') ;    %  then change to input value
    iforcenaivesp = getfield(paramstruct,'iforcenaivesp') ; 
  end ;

  if isfield(paramstruct,'icolor') ;    %  then change to input value
    icolor = getfield(paramstruct,'icolor') ; 
  end ;

  if isfield(paramstruct,'markerstr') ;    %  then change to input value
    markerstr = getfield(paramstruct,'markerstr') ; 
  end ;

  if isfield(paramstruct,'isubpopkde') ;    %  then change to input value
    isubpopkde = getfield(paramstruct,'isubpopkde') ; 
  end ;

  if isfield(paramstruct,'idataconn') ;    %  then change to input value
    idataconn = getfield(paramstruct,'idataconn') ; 
  end ;

  if isfield(paramstruct,'idataconncolor') ;    %  then change to input value
    idataconncolor = getfield(paramstruct,'idataconncolor') ; 
  end ;

  if isfield(paramstruct,'idataconntype') ;    %  then change to input value
    idataconntype = getfield(paramstruct,'idataconntype') ; 
  end ;

  if isfield(paramstruct,'ibigdot') ;    %  then change to input value
    ibigdot = getfield(paramstruct,'ibigdot') ; 
  end ;

  if isfield(paramstruct,'idatovlay') ;    %  then change to input value
    idatovlay = getfield(paramstruct,'idatovlay') ; 
  end ;

  if isfield(paramstruct,'ndatovlay') ;    %  then change to input value
    ndatovlay = getfield(paramstruct,'ndatovlay') ; 
  end ;

  if isfield(paramstruct,'datovlaymax') ;    %  then change to input value
    datovlaymax = getfield(paramstruct,'datovlaymax') ; 
  end ;

  if isfield(paramstruct,'datovlaymin') ;    %  then change to input value
    datovlaymin = getfield(paramstruct,'datovlaymin') ; 
  end ;

  if isfield(paramstruct,'legendcellstr') ;    %  then change to input value
    legendcellstr = getfield(paramstruct,'legendcellstr') ; 
  end ;

  if isfield(paramstruct,'mlegendcolor') ;    %  then change to input value
    mlegendcolor = getfield(paramstruct,'mlegendcolor') ; 
  end ;

  if isfield(paramstruct,'maxlim') ;    %  then change to input value
    maxlim = getfield(paramstruct,'maxlim') ; 
  end ;

  if isfield(paramstruct,'iplotaxes') ;    %  then change to input value
    iplotaxes = getfield(paramstruct,'iplotaxes') ; 
  end ;

  if isfield(paramstruct,'iplotdirvec') ;    %  then change to input value
    iplotdirvec = getfield(paramstruct,'iplotdirvec') ; 
  end ;

  if isfield(paramstruct,'ibelowdiag') ;    %  then change to input value
    ibelowdiag = getfield(paramstruct,'ibelowdiag') ; 
  end ;

  if isfield(paramstruct,'titlecellstr') ;    %  then change to input value
    titlecellstr = getfield(paramstruct,'titlecellstr') ; 
  end ;

  if isfield(paramstruct,'titlefontsize') ;    %  then change to input value
    titlefontsize = getfield(paramstruct,'titlefontsize') ; 
  end ;

  if isfield(paramstruct,'labelcellstr') ;    %  then change to input value
    labelcellstr = getfield(paramstruct,'labelcellstr') ; 
  end ;

  if isfield(paramstruct,'labelfontsize') ;    %  then change to input value
    labelfontsize = getfield(paramstruct,'labelfontsize') ; 
  end ;
  
  if isfield(paramstruct,'savestr') ;    %  then use input value
    savestr = getfield(paramstruct,'savestr') ; 
    if  ~ischar(savestr)  &  ~isempty(savestr) ;   
                          %  then invalid input, so give warning
      disp('!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!') ;
      disp('!!!   Warning from scatplotSM.m:   !!!') ;
      disp('!!!   Invalid savestr,             !!!') ;
      disp('!!!   using default of no save     !!!') ;
      disp('!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!') ;
      savestr = [] ;
    end ;
  end ;

  if isfield(paramstruct,'iscreenwrite') ;    %  then change to input value
    iscreenwrite = getfield(paramstruct,'iscreenwrite') ; 
  end ;

  if isfield(paramstruct,'celltypes' );  % then change to input value
    celltypes = paramstruct.celltypes;
  end ;
  
  if isfield(paramstruct,'cellsubtypes' );  % then change to input value
    cellsubtypes = paramstruct.cellsubtypes;
  end ;
  
  if isfield(paramstruct,'slidenames' );  % then change to input value
    slidenames = paramstruct.slidenames;
  end ;
  
  if isfield(paramstruct,'slideids' );  % then change to input value
    slideids = paramstruct.slideids;
  end ;



end ;    %  of resetting of input parameters

%  Set up output for mexplorer
%  i.e.  set flag to augment with user data
%
if ( ~isempty( slidenames ) & ~isempty( slideids ) & ~isempty( celltypes ) )
  augmentWithUserData = 1;
  if ( isempty( cellsubtypes ) )
    cellsubtypes = cell( size( celltypes ) ); % just initialize it empty
  end
else
  augmentWithUserData = 0;
end

%  set preliminary stuff
%
d = size(data,1) ;
         %  dimension of each data curve
n = size(data,2) ;
         %  number of data curves
if ~isempty(mdir) ;
  ncomp = size(mdir,2) ;
           %  number of components to plot
  if ~(d == size(mdir,1)) ;
    disp('!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!') ;
    disp('!!!   Error from scatplotSM.m:    !!!') ;
    disp('!!!   Dimension mdir must be      !!!') ;
    disp('!!!   same as dimension of data   !!!') ;
    disp('!!!   Terminating execution       !!!') ;
    disp('!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!') ;
    return ;
  end ;
else ;    %  no direction vectors entered, so are doing pca
  ncomp = 0 ;
  if npcadiradd < 1 ;
    disp('!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!') ;
    disp('!!!   Warning from scatplotSM.m:                             !!!') ;
    disp('!!!   mdir is empty, so will assume PCAdirections desired    !!!') ;
    disp('!!!   but npcadiradd needs to be > 0                         !!!') ;
    disp('!!!   Will reset to 4 PC directions                          !!!') ;
    disp('!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!') ;
    npcadiradd = 4 ;
  end ;
end ;

if (ncomp + abs(npcadiradd)) < 2 ;
  disp('!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!') ;
  disp('!!!   Error from scatplotSM.m:                               !!!') ;
  disp('!!!   # Columns of mdir + abs(npcadiradd) must be >= 2       !!!') ;
  disp('!!!   For 1-d Projections, use projplot1dSM.m                !!!') ;
  disp('!!!   Terminating execution                                  !!!') ;
  disp('!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!') ;
  return ;
end ;

if ~isempty(maxlim) ;
  if ~(maxlim == 1) ;
    if  ~(size(maxlim,1) == (size(mdir,2) + abs(npcadiradd)))  | ...
        ~(size(maxlim,2) == 2)  ;
      disp('!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!') ;
      disp('!!!   Warning from scatplotSM.m:                   !!!') ;
      disp(['!!!   Invalid size of maxlim: ' num2str(size(maxlim,1)) ' x ' num2str(size(maxlim,2))]) ;
      disp(['!!!   Must be: ' num2str(size(mdir,2) + abs(npcadiradd)) ' x 2']) ;
      disp('!!!   Resetting maxlim to default of []            !!!') ;
      disp('!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!') ;
      maxlim = [] ;
    end ;
  end ;
end ;

if  (size(icolor,1) > 1)  |  (size(icolor,2) > 1)  ;    %  if have color matrix
  if ~(3 == size(icolor,2)) ;
    disp('!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!') ;
    disp('!!!   Warning from scatplotSM.m:               !!!') ;
    disp('!!!   icolor as a matrix must have 3 columns    !!!') ;
    disp('!!!   Resetting to icolor = 1, Matlab default   !!!') ;
    disp('!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!') ;
    icolor = 1 ;
  elseif ~(n == size(icolor,1)) ;
    disp('!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!') ;
    disp('!!!   Warning from scatplotSM.m:               !!!') ;
    disp(['!!!   icolor as a matrix must have ' num2str(n) ' rows']) ;
    disp('!!!   Resetting to icolor = 1, Matlab default   !!!') ;
    disp('!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!') ;
    icolor = 1 ;
  end ;
end ;

icolorprint = 1 ;
if  size(icolor,1) == 1  &  size(icolor,2) == 1 ;
  if isstr(icolor) ;
    if strcmp(icolor,'k') ;
      icolorprint = 0 ;
    end ;
  elseif icolor == 0  ;
    icolorprint = 0 ;
  end ;
end ;



%  Setup labelling
%
if isempty(labelcellstr) ;
  vlabelflag = ones(ncomp,1) ;
      %  use default labelling for all
else ;
  vlabelflag = [] ;
  for ilabel = 1:length(labelcellstr) ;
    if isempty(labelcellstr(ilabel)) ;
      vlabelflag = [vlabelflag; 0] ;
          %  use no label in this direction
    else ;
      vlabelflag = [vlabelflag; 2] ;
          %  use given label in this direction
    end ;
  end ;

  if length(labelcellstr) < ncomp ;
    vlabelflag = [vlabelflag; ones(ncomp - length(labelcellstr),1)] ;
        %  fill in with default labels for the rest
  end ;

end ;



%  Mean Recenter Data (if needed)
%
if irecenter == 0 ;    %  Do not do mean recentering
  datac = data ;
  if iscreenwrite == 1 ;
    disp('    scatplotSM:    Skipping Mean Recentering') ;
  end ;
else ;    %  Mean Recenter Data
  datac = data - vec2matSM(mean(data,2),n) ;
  if iscreenwrite == 1 ;
    disp('    scatplotSM:    Finished Mean Recentering') ;
  end ;
end ;



%  Add some PCA directions (if needed)
%
if abs(npcadiradd) > 0 ;    %  then add PCA Direction vectors to mdir

  if iscreenwrite == 1 ;
    disp('    scatplotSM:    Adding PC Directions') ;
  end ;

  if npcadiradd > 0 ;    %  then add full data PC dir'ns
    datasub = datac ;
  else ;    %  add PC dirn's in subspace orthogonal to mdir
    %    mproj = mdir * pinv(mdir' * mdir) * mdir' ;
    %    datasub = data - mproj * data ;
          %  The above were old lines, which worked,
          %  but would explode memory in really high
          %  dimensional situations.  
          %  Problem was that  mdir' * mdir  created too large a matrix
          %  Avoid this with a careful SVD calculation:
    r = rank(mdir) ; 
    [U,S,V] = svd(mdir,'econ') ;
    U = U(:,1:r) ;
        %  orthonormal basis matrix of mdir column space
    mprojdata = U * (U' * datac) ;
        %  data projected onto mdir column space
    datasub = datac - mprojdata ;
        %  data projected onto orthogonal space
  end ;
  PCAparamstruct = struct('npc',abs(npcadiradd), ...
                          'iorient',iorient, ...
                          'viout',[0 1]) ;
        % iorient = 3 tends to put "scores in right quadrant"
  outstruct = pcaSM(datasub,PCAparamstruct) ;
  meigvec = getfield(outstruct,'meigvec') ;


  npcada = size(meigvec,2) ;
  if npcada < abs(npcadiradd) ;
    disp('!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!') ;
    disp('!!!   Warning from scatplotSM.m:   !!!') ;
    disp('!!!   Added PCs not of full rank   !!!') ;
    disp(['!!!   Reducing to ' num2str(npcada) ' added PCs']) ;
    disp('!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!') ;
  end ;

  mdir = [mdir meigvec] ;

  if length(labelcellstr) > ncomp ;    %  then there are input labels to use

    for ilabel = (ncomp + 1):length(labelcellstr) ;
      if isempty(labelcellstr{ilabel}) 
        vlabelflag = [vlabelflag; 0] ;
            %  use no label in this direction
      else ;
        vlabelflag = [vlabelflag; 2] ;
            %  use given label in this direction
      end ;
    end ;

    if length(labelcellstr) < ncomp + abs(npcadiradd) ;
      vlabelflag = [vlabelflag; 2 * ones(ncomp - length(labelcellstr),1)] ;
          %  fill in with "pc" labels for the rest

      if length(labelcellstr) < ncomp ;
          %  then need to pad labelcellstr with {}s
        ilabelstart = length(labelcellstr) + 1 ;
        for ilabel = ilabelstart:ncomp ;
          labelcellstr = cat(1,labelcellstr,{}) ;
        end ;
      end ;

      for ipc = 1:abs(npcadiradd) ;
        if npcadiradd > 0 ;
          labelcellstr = cat(1,labelcellstr,{['PC' num2str(ipc) ' Direction']}) ;
        else ;
          labelcellstr = cat(1,labelcellstr,{['Ortho PC' num2str(ipc) ' Direction']}) ;
        end ;
      end ;

    end ;

  else ;    % use default labels

    vlabelflag = [vlabelflag; 2 * ones(abs(npcadiradd),1)] ;
    sizelcs = size(labelcellstr,1) ;
    if sizelcs < ncomp ;    %  then need to pad with empty elements
      for i = 1:(ncomp - sizelcs) ;
        labelcellstr = cat(1,labelcellstr,{''}) ;
      end ;
    end ;
    for ipc = 1:abs(npcadiradd) ;    %    fill in with "pc" labels for the rest
      if npcadiradd > 0 ;
        labelcellstr = cat(1,labelcellstr,{['PC' num2str(ipc) ' Direction']}) ;
      else ;
        labelcellstr = cat(1,labelcellstr,{['Ortho PC' num2str(ipc) ' Direction']}) ;
      end ;
    end ;

  end ;


  ncomp = ncomp + abs(npcada) ;
         %  number of components to plot


end ;



%  make draftsmans' 2-d plot
%
if iscreenwrite == 1 ;
  disp('    scatplotSM:    Making 2d Draftsman''s Plots') ;
end ;
clf ;


%  First work on diagonal, and put down same kde stuff
%
for ic = 1:ncomp ;

  if size(titlecellstr,2) >= ic ;
    if size(titlecellstr{ic},1) >= ic ;
      titlestr = titlecellstr{ic}(ic,:) ;
    else ;
      titlestr = '' ;    
    end ;
  else ;
    titlestr = '' ;    
  end ;

  if vlabelflag(ic) == 0 ;
    xlabstr = '' ;
  elseif vlabelflag(ic) == 1 ;
    xlabstr = ['Direction ' num2str(ic)]  ;
  else ;
    xlabstr = labelcellstr{ic} ;
  end ;

  if isempty(maxlim) ;
    vaxlim1 = [] ;
  elseif maxlim == 1 ;
    vaxlim1 = 1 ;
  else ;
    vaxlim1 = maxlim(ic,:) ;
  end ;

  subplot(ncomp,ncomp,(ncomp+1)*(ic-1)+1) ;

    paramstruct1 = struct('icolor',icolor, ...
                          'markerstr',markerstr, ...
                          'isubpopkde',isubpopkde, ...
                          'ibigdot',ibigdot, ...
                          'idatovlay',idatovlay, ...
                          'ndatovlay',ndatovlay, ...
                          'datovlaymax',datovlaymax, ...
                          'datovlaymin',datovlaymin, ...
                          'vaxlim',vaxlim1, ...
                          'titlestr',titlestr, ...
                          'titlefontsize',titlefontsize, ...
                          'xlabelstr',xlabstr, ...
                          'labelfontsize',labelfontsize, ...
                          'ifigure',0, ...
                          'savestr',[], ...
                          'iscreenwrite',iscreenwrite) ;
    if ic == 1 ;
      if ~(isempty(legendcellstr)) ;
        paramstruct1 = setfield(paramstruct1,'legendcellstr',legendcellstr) ;
      end ;
      if ~(isempty(mlegendcolor)) ;
        paramstruct1 = setfield(paramstruct1,'mlegendcolor',mlegendcolor) ;
      end ;
    end ;

    if ( augmentWithUserData )
      paramstruct1.celltypes = celltypes;
      paramstruct1.cellsubtypes = cellsubtypes;
      paramstruct1.slidenames = slidenames;
      paramstruct1.slideids = slideids;
    end

    projplot1SM(datac,mdir(:,ic),paramstruct1) ;


end ;    %  of loop for diagonals



%  Now do off diagonals
%
for ic = 1:ncomp ;

  if vlabelflag(ic) == 0 ;
    ylabstr = '' ;
  elseif vlabelflag(ic) == 1 ;
    ylabstr = ['Direction ' num2str(ic)]  ;
  else ;
    ylabstr = labelcellstr{ic} ;
  end ;


  if ibelowdiag == 0 ;
    jstart = ic + 1 ;
  else ;
    jstart = 1 ;
  end ;
  for jc = jstart:ncomp ;

    if ic ~= jc ;   %  then not on diagonal, so make a plot

      if size(titlecellstr,2) >= jc ;
        if size(titlecellstr{jc},1) >= ic ;
          titlestr = titlecellstr{jc}(ic,:) ;
        else ;
          titlestr = '' ;
        end ;
      else ;
        titlestr = '' ;
      end ;

      if vlabelflag(jc) == 0 ;
        xlabstr = '' ;
      elseif vlabelflag(jc) == 1 ;
        xlabstr = ['Direction ' num2str(jc)]  ;
      else ;
        xlabstr = labelcellstr{jc} ;
      end ;

      if isempty(maxlim) ;
        vaxlim2 = [] ;
      elseif maxlim == 1 ;
        vaxlim2 = 1 ;
      else ;
        vaxlim2 = [maxlim(jc,:) maxlim(ic,:)] ;
      end ;


      subplot(ncomp,ncomp,ncomp*(ic-1)+jc) ;
      
        paramstruct2 = struct('icolor',icolor, ...
                              'markerstr',markerstr, ...
                              'iforcenaivesp',iforcenaivesp, ...
                              'idataconn',idataconn, ...
                              'idataconncolor',idataconncolor, ...
                              'idataconntype',idataconntype, ...
                              'ibigdot',ibigdot, ...
                              'vaxlim',vaxlim2, ...
                              'iplotaxes',iplotaxes, ...
                              'iplotdirvec',iplotdirvec, ...
                              'titlestr',titlestr, ...
                              'titlefontsize',titlefontsize, ...
                              'xlabelstr',xlabstr, ...
                              'ylabelstr',ylabstr, ...
                              'labelfontsize',labelfontsize, ...
                              'ifigure',0, ...
                              'savestr',[], ...
                              'iscreenwrite',iscreenwrite) ;

        if ( augmentWithUserData )
          paramstruct2.celltypes = celltypes;
          paramstrcut2.cellsubtypes = cellsubtypes;
          paramstruct2.slidenames = slidenames;
          paramstruct2.slideids = slideids;
        end
        
        projplot2SM(datac,[mdir(:,jc) mdir(:,ic)],paramstruct2) ;


    end ;

  end ;

end ;




if ~isempty(savestr) ;   %  then create postscript file

  orient landscape ;

  if icolorprint ~= 0 ;     %  then make color postscript
    print('-dpsc',savestr) ;
  else ;                %  then make black and white
    print('-dps',savestr) ;
  end ;
  
  if ( augmentWithUserData )
    saveas(gcf, savestr, 'fig')
  end ;

end ;


