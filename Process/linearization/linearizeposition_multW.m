function [statematrix, segmenttable, trajwells, wellSegmentInfo, segmentInfo] = linearizeposition_multW(directoryname,fileprefix, index, varargin)

% Originally Mattias Karlsson, edited Kenny Kay, large sections modified for
% multi-W track by Mari Sosa, January 2017. 

% Takes each position from the animal's 'pos' structure  and 
% projects it onto the linear trajectory segments for the environment
% and task defined by coordprogram.
% This routine finds the segment closest to each position, projects
% the position to the closest location on the segment. It then calculates
% the linear distance from each well to the animal and which well-to-well
% trajectory the animal is on for each time step.
%
%      --INPUTS--
%
%    directoryname - example 'data99/user/animaldatafolder/', a folder 
%                    containing processed matlab data for the animal
%       fileprefix - animal specific prefix for each datafile
%            index - [day epoch]
%
%          options - 'lowercasethree' - varaible prefix (default '')
%
%                    'maxvelocity' - maximum allowed velocity in cm/sec (default ~300 for W-track, 500 for multi-W)
%                                  - can input as a varargin to linearizepos
%
%                    'welldist' - to decide whether or not the animal
%                    completed the trajectory, the program needs a radius
%                    around the endpoint of each trajectory which the
%                    animal must enter to trigger a completed trajectory (default 5 cm for W-track)
%
%                    'mindiff' - the minimum amount of time between two
%                    foodwell zone entries to be concidered a real
%                    trajectory (default 2 seconds)
%
%                    'velocitysmoothwindow' - to compute linear speed, the
%                    program needs to smooth linear position data.  It is smoothed
%                    with a gaussian of length VSW and std VSW/4 (default 2
%                    seconds)
%
%                    'welldistthresh' -the threshold linear distance (in cm) from the well for the program 
%                    to add a trajectory when the animal turns around
%                    before reaching the end well.  (default 0)
%
%                   'branchpointdist' -this performs a similar function to
%                   welldistthresh, except the distance is defined as the
%                   distance traveled a the part of the track that only leads to one well.
%                   If a value is given for this, it supercedes any value
%                   given for welldistthresh. If set to Inf, welldisthresh
%                   supercedes.
%          
%                   'manual_segment' - [day1 epoch1 startindex1 endindex1 segmentnumber1 ;  day2 epoch2 ...]
                        %   each row corresponds to a set of indices that
                        %   need to be manually changed -- this is after
                        %   you've verified separately that these indices
                        %   are currently being misclassified
                        %   (dfskk_thetaprecess)
                        % added kk 3.24.14, not implemented
%      
%
%       --OUTPUTS--
%
%      statematrix - a structure with each field the same length as pos.data
%                    Contains information about the linear position of the animal
%                    for each time step. 
%
%     segmenttable - a lookup table with 3 columns: the first column
%                    is the linear trajectory number (the row number in 
%                    lindistpos{i}{e}.traject), the second column is the segment number in 
%                    that trajectory, and the third column is the unique indentifier number 
%                    for that segment, which is used in the segmentIndex field of statematrix. Any segments 
%                    that are used in multiple trajectories are given the same identifier number.
%                     
%        trajwells - gives the start and end well numbers for each trajectory (each row is 
%                    a trajectory).  Each trajecory is defined from the
%                    output of coordprogram.
%
%  wellSegmentInfo - a structure with fields describing aspects of the
%                    segments containing a reward well, such as the well's distance to an
%                    intersection and each segment's direction relative to the wells
%
%      segmentInfo - a structure with fields describing each segment, such
%                    as segment length and which segments are connected to each other


manual_segment = [];
maxv = 500;
lowercasethree = ''; %default variable prefix is none
welldist = 9;
mindiff = 2;
welldistthresh = 0;
branchpointdist = Inf;
smoothwidth = 2;
if (length(varargin) == 1)
    varargin = varargin{1};
end
%set variable options
for option = 1:2:length(varargin)-1
    
    switch varargin{option}
        case 'lowercasethree'
            lowercasethree = varargin{option+1};
        case 'maxvelocity'
            maxv = varargin{option+1};          
        case 'maxsegdiff'
            maxsegdiff = varargin{option+1};          
        case 'welldist'
            welldist = varargin{option+1};
        case 'mindiff'
            mindiff = varargin{option+1};
        case 'velocitysmoothwindow'
            smoothwidth = varargin{option+1};
        case 'welldistthresh'
            welldistthresh = varargin{option+1};
        case 'branchpointdist'
            branchpointdist = varargin{option+1};
        case 'manual_segment'
            manual_segment = varargin{option+1};       % 3.24.14 kk added
        otherwise
            error(['Option ''', varargin{option}, ''' not defined']);
    end
end
    
dsz = '';
if (index(1) < 10)
   dsz = '0';
end

%load the data
eval(['load ',directoryname,fileprefix,'pos', dsz, num2str(index(1)), '.mat']);
eval(['pos = ',lowercasethree,'pos;'])
eval(['load ',directoryname,fileprefix,'task',dsz, num2str(index(1)), '.mat']);
eval(['task = ',lowercasethree,'task;'])


% make sure we have the direction information
toknum = isdatafield(pos{index(1)}{index(2)}.fields, 'dir-sm');
if (~toknum)
    disp('No dir-sm direction field in pos! try obtaining it -- continuing w/ ordinary dir field');
    toknum = isdatafield(pos{index(1)}{index(2)}.fields, 'dir');
    if ~toknum
        error('No dir-sm direction field in pos!')
    end
end

pos = pos{index(1)}{index(2)}.data;
task = task{index(1)}{index(2)};
timestep = pos(2,1) - pos(1,1);

%initialize variables
poslength = size(pos,1);
segment = zeros(poslength,1);
lindist = ones(size(pos,1), 1) * -1;
vect = zeros(poslength,2);
newpos = pos;

coordInTime = task.linearcoord;
% clear task;
% if no taskstruct, use the coordinate program to fetch the createtaskstruct (getcoord) track coordinates 
%eval(['coordProgramHandle = @',coordprogram,';']);
%coordInTime = feval(coordProgramHandle,pos(:,1),task);

%distsum: calculate the total linear distance to each coordinate
      % {i} is trajectory
      % each row is a distance to one of the clicked-on coordinates 
distsum = [];
for i = 1:length(coordInTime)
    firstcoord{i} = coordInTime{i}(:,:,1);
    if (size(firstcoord{i},1) > 1)
        distsum{i}(1,1) = 0;
        for j = 2:size(firstcoord{i},1)
            distsum{i}(j,1) = sqrt( ((firstcoord{i}(j,1) - firstcoord{i}(j-1,1))^2) + ((firstcoord{i}(j,2) - firstcoord{i}(j-1,2))^2) ) + ...
                distsum{i}(j-1);
        end
    end
end

%create a table of segments describing which trajectory and segnum
%they are.  This table give segment an identifying number. If any segment
%is used in multiple trajectories, it is only listed once as belonging to
%the first trajectory. This function also calculates specific information
%about the track segments and trajectories. 
[segmenttable, wellindex, trajwells, wellSegmentInfo, segmentInfo] = getSegmentTable(firstcoord);

    %segmenttable :
        % 1st col: trajectory # (1 or 2) - defined in getcoord_<x>track
        % 2nd col: trajectory segment # (1st, 2nd, or 3rd segment OF THAT trajectory)
        % 3rd col: unique segment # (1,2,3,4,5 -- for W-track at least)
    % wellindex:
        % (one row for each well)
        % (note that there is no need to have 4 rows.. since [1 1] suffices for the home well --> for W track)
        % 1st col is trajectory # (1 or 2)
        % 2nd col is coordinate number (1, 2, 3, or 4 -- the 4 coords of a trajectory)
    %trajwells: pairs of wells that are connected on each of the trajectories
    %wellSegmentInfo:
        % .distanceTable: distance from the start of each segment to each
            % of the wells (columns are seg numbers, if there is a 0 in the
            % column, then that segment starts with a well; 
        % .segmentIndex: unique segment #s that are well arms; in mult W, segs 1,3,5,7,9,11
        % .distantToIntersection: length of the well arms
        % .wellCoord = well coordinates as clicked in getcoord
        % .pathTable = table of every possible path between n connected segments
        % .segmentDirection = table of which direction the segment will
            % linearize in depending on its starting point (1 forward, 0
            % backward)
    %segmentInfo:
        % .segmentCoords = [x y x y] of the 2 points that bookend each segment, each row is a segment
        % .segmentLength = length of each segment in centimeters
        % .connectedSegments = <traj#> cells, in which each cell contains the segment IDs directly connected to the segment <traj#>
        
ntraj = size(trajwells,1);

%get the well locations for all time frames
for i = 1:size(wellindex,1)
    welllocations(i,1:2,:) = coordInTime{wellindex(i,1)}(wellindex(i,2),1:2,:);
end

for i = 1:size(pos,1)
    % iterate through each unique segment # (1,2,3,4,5, etc) 
    for findcoord = 1:length(segmentInfo.segmentLength)
        % tableInd is the row in segmenttable corresponding to the unique segment # 
        tableInd = min(find(segmenttable(:,3) == findcoord));  
            % here, obtain the coordinates of the START (1st row) and END (2nd row) points of each segment
        coord{findcoord} = coordInTime{segmenttable(tableInd,1)}(segmenttable(tableInd,2):segmenttable(tableInd,2)+1,:,i);
        %calculate each segment's vector
        coordvector(findcoord,1:2) = diff(coord{findcoord});
            % coordvector:
                % each row corresponds to each unique segment # (1,2,3,4,5)
                % [xend-xstart yend-ystart]  
    end
    
end

% in linearizepos the main loop projects each 2D position point onto each segment (and does error correction along the way).
% will call supporting function projectpoint
[linpos, segmentIndex, segdist] = createlinpos(pos, task, index, coord);

%  fill in any missing elements skipped from invalid positions
segmentIndex = vectorfill(segmentIndex,0);

% at all time points, find the linear distance from each well
for s = 1:length(segmentIndex)
    newsegment = segmentIndex(s);
    % find the vector along each segment
    vect(s,1) = coordvector(newsegment,1);  % x distance between the current trajectory's begin x and end x
    vect(s,2) = coordvector(newsegment,2);  % y distance between the current trajectory's begin y and end y  
            %find the linear distance for the point from each well
            %this requires a check for which direction the segment is
            %aligned relative to the well
            distToSegment = wellSegmentInfo.distanceTable(:,newsegment)';  % the distance from EACH well (1,2,3, in columns) to the newSegment in question
            segmentDist = [];
            % iterate over wells
            % (wellSegmentInfo.distanceTable -- is wellno x segment #)
                    % installs a 1 if the well is POTENTIALLY at the endpoint of the segment, otherwise 0
            for wellcount = 1:size(wellSegmentInfo.distanceTable,1)
                if  (wellSegmentInfo.segmentDirection(wellcount,newsegment) == 1)
                    % == 1 : it is aligned in the forward direction, then directly report segdist
                    segmentdist(1,wellcount)  = segdist(s);
                else
                    % == 0 : it is aligned in the backward direction (?), then report (segment length - segdist)
                    segmentdist(1,wellcount)  = segmentInfo.segmentLength(newsegment)-segdist(s);
                end
            end
            % segmentdist (three values, one for each well):
                % each entry is distance along the segment -- with well target taken into account
                    % if a well is being traveled towards is the "segdist"
                    % if a well is being traveled away from, then it's
                    % segment length - segdist
            
            %calculate the linear distance from each well
            lindist(s,1:wellcount) = distToSegment + segmentdist;   
end

%  fill in any missing elements skipped from invalid positions
segdist = vectorfill(segdist,-1);
vect(:,1) = vectorfill(vect(:,1),0);
vect(:,2) = vectorfill(vect(:,2),0);

for i = 1:size(lindist,2)
    lindist(:,i) = vectorfill(lindist(:,i),-1);
end

%get the foodwell enter and exit times: trajmatrix = [starttime startwell endwell]; starttime is exit time from start well
trajmatrix = gettraject(pos, welllocations,trajwells, lindist, wellSegmentInfo, welldist, mindiff, welldistthresh,branchpointdist);

%%%%%%% Now, we want to find out which well to use as the reference point for linear distance traveled. %%%%%%%
% This depends on which trajectory the animal is on.
whichTraj = zeros(size(pos,1),1);

% (Not for Multi-W) if all trajectories start with the same well, then we should always use
%that well as the reference point.
if ( (sum(abs(diff(trajwells(:,1)))) == 0) | (isempty(diff(trajwells(:,1)))) )
    whichTraj(:) = trajwells(1,1); %just pick the first traj for all points, later the reference well will be chosen as the trajectory's starting point
end
if ~isempty(trajmatrix)
        %for each timestamp in pos, get which of these complete trajectories the
        %animal is in
        welltraj = trajmatrix(lookup(pos(:,1),trajmatrix(:,1),-1),2:3);
else
        disp('WARNING: No completed trajectories.  Check endpoint locations if this is incorrect.')
        %because no trajectories were competed, we fill the start and end wells
        %with NaN's
        welltraj = nan(size(pos,1),2);
end
samewelltraj = [];

if (whichTraj(1) == 0) %more than one starting well exists, so we need more complicated work to find the reference well
    
    % (Not for Multi-W) if the animal was on a segment that belongs to only one trajectory,
    %then we assign that trajectory to that time
    singleTrajSegments = find(rowcount(segmentIndex,segmenttable(:,3))==1);
    whichTraj(singleTrajSegments) = segmenttable(rowfind(segmentIndex(singleTrajSegments),segmenttable(:,3)),1);
    zeroind = find(whichTraj == 0); %these are the indeces to the times when the animal was on amiguous segments

    % (FOR MULTI-W) next, we look at the times when the animal was on segments belonging
    %to 2 or more trajectories.  We use the start and end wells for each
    %time step to determine the best trajectory to assign the timestep to.
    
    if ~isempty(trajmatrix)
    badindex=[]; % initialize list of indices that could not be assigned a complete trajectory  
    
        %To find the reference well (below), we will temporarily find the trajectory number for each well exit-to-entry, if it exists.
        %To do this we must check both the forward and reverse directions.
        %This is not the final trajectory assignment, but 'bad' indices
        %will be accurate.
        tmptrajnum = rowfind(welltraj(zeroind,:),trajwells);
        tmptrajnum2 = rowfind(welltraj(zeroind,:),trajwells(:,[2 1])); %zeroind
        tmptrajnum = max(tmptrajnum,tmptrajnum2);
        trajNonzeros = find(tmptrajnum > 0);
        %only include the times when the segment belonged to the trajectory
        segNonzero = segmentIndex(trajNonzeros);
        nonstandard = zeros(size(whichTraj)); % any flag >0 will indicate nonstandard, and will indicate the number of segments deviated
        
        % Loop through the times with identified possible trajectories and
        % check to see if the current segment is a member of that traj
        i=1;
        while i<=length(segNonzero)
            segmentsInTraj = segmenttable(segmenttable(:,1)==tmptrajnum(trajNonzeros(i)),3);
            if ismember(segNonzero(i),segmentsInTraj) % if the segment is a member, assign this time point to the trajectory
                whichTraj(trajNonzeros(i)) = tmptrajnum(trajNonzeros(i));
                i = i+1;
            else %if current segment is not a member of the trajectory
                % check if this is an excursion that still results in the completion of the trajectory by checking neighboring segments 
                if i>1 && ismember(segNonzero(i),segmentInfo.connectedSegments{segNonzero(i-1)})
                    currsegment = segNonzero(i);
                    new = segNonzero(i:end);
                    % find the next index that contains a segment included in the trajectory, and is still tagged as the same trajectory
                    nextsegind = find(ismember(new,segmentsInTraj) & tmptrajnum(trajNonzeros(i:end))==tmptrajnum(trajNonzeros(i)),1);
                    % if all timesteps between now and the next seg that satisfied the above are the same as the current seg, we transcribe the trajectory
                    if (~isempty(nextsegind) && (sum((new(1:nextsegind-1) == currsegment) == 0) == 0))
                        inds = (1:nextsegind-1)+(i-1); %convert the "new" indices back to the ones we've been using
                        whichTraj(trajNonzeros(inds)) = tmptrajnum(trajNonzeros(inds));
                        nonstandard(trajNonzeros(inds)) = 1; %deviation by 1 segment
                        i = inds(end)+1;
                    % if a nextseg was found but there is more than one segment that is not a member of the current trajectory, still transcribe    
                    elseif (~isempty(nextsegind) && (sum((new(1:nextsegind-1) == currsegment) == 0) > 0))
                        inds = (1:nextsegind-1)+(i-1); %convert the "new" indices back to the ones we've been using
                        whichTraj(trajNonzeros(inds)) = tmptrajnum(trajNonzeros(inds));
                        numsegs = sum(unique(new(1:nextsegind-1))~=currsegment);
                        nonstandard(trajNonzeros(inds)) = numsegs; %deviation by numsegs
                        i = inds(end)+1;
                    % if we could not find another segment that's part of this trajectory before the end of the data, we skip the index
                    % this is probably because the trajectory was not completed
                    else
                        badindex = [badindex; i];
                        i = i+1; % continue without replacing values in whichTraj
                    end
                % if the current segment is not even connected to the
                % segments of this trajectory, skip the index.
                % NOTE: if there is maze exploration before the first well
                % entry, those indices will be unassigned!!
                else
                    badindex = [badindex; i];
                    i = i+1; % continue without replacing values in whichTraj
                end
            end
        end
        disp(sprintf('%d indices were not assigned a trajectory',length(badindex)))

        zeroind = find(whichTraj == 0);
        % find trajectories where start and end well are the same (i.e. trackback)
        samewelltraj = welltraj(zeroind(welltraj(zeroind,1) == welltraj(zeroind,2)),1);  
        samewelltrajind = zeroind(welltraj(zeroind,1) == welltraj(zeroind,2));
    end   
end

referencewell = ones(size(pos,1),1) *-1;
nonzeros = find(whichTraj>0);
% assign reference well as the start well of the trajectory
referencewell(nonzeros) = trajwells(whichTraj(nonzeros),1);
%assign the times when the animal exited and entered the same well to that
%reference well (if the time step has not already been defined to a
%reference well) (i.e. trackbacks)
if ~isempty(samewelltraj)
    referencewell(samewelltrajind) = samewelltraj;
end
%any remaining -1's in the referencewell vector are times that are
%too difficult to assign a reference well to 

%calculate the segment direction for each time step
segmentdir = ones(size(pos,1),size(lindist,2));
    % iterate through each well
for i = 1:size(lindist,2)
    % for each well (col), determine at each time point (row) whether the
    % animal is moving toward (1) or away from well
    segmentdir(:,i) = wellSegmentInfo.segmentDirection(i,segmentIndex);
end


% compute head directions relative to the track segments (positive is moving toward the right of the maze (as defined by click order in getcoord_6armtrack),
%  negative is moving toward the left of the maze (reverse direction of each traj defined by getcoord) 
vlen = sqrt(vect(:,1).^2 + vect(:,2).^2);   % vector lengths of each segment's START pt to END pt vector
% normalize vect
vect(:,1) = vect(:,1) ./ vlen;              % at this point, vect is the (x,y) list of normal vectors pointing from the START pt to END pt of each current segment
vect(:,2) = vect(:,2) ./ vlen;
headdir = pos(:,toknum);
y = sin(headdir);                           % this x, y vector is the set of normal vectors of the animal's actual current head direction
x = cos(headdir);
segheaddir = (vect(:,1) .* x + vect(:,2) .* y);      % if pointing in same direction, = 1  // otherwise = -1
% simply copy the vector over for each of the other wells
segheaddir = repmat(segheaddir,[1 size(lindist,2)]);
%reverse the head direction for the times when the segment linearization is
%reverse relative to the trajectory
reverseinds = find(segmentdir == 0);
segheaddir(reverseinds) = segheaddir(reverseinds)*-1;


%Compute velocity
npoints = smoothwidth/timestep;
filtstd = smoothwidth/(4*timestep);
% the default filter for smoothing motion direction is a n second long gaussian
% with a n/4 second stdev
filt = gaussian(filtstd, npoints);
% smooth the linear distances with the filter and then go through the linear
% distance positions and take all of the differences to get velocities
velocity = [];
for i = 1:size(lindist,2)
    smoothdist = smoothvect(lindist(:,i), filt);
    v = diff(smoothdist) / timestep;
    v = [v(1); v];
    velocity = [velocity v];
end

statematrix.time = pos(:,1);                    % time
statematrix.segmentIndex = segmentIndex;        % segment number rat is on at each timestep
statematrix.wellExitEnter = welltraj;           % start and end well that bookend the trajectory the rat is on at each timestep
statematrix.segmentHeadDirection = segheaddir;  % rat's head direction relative to the direction of the segment
statematrix.linearVelocity = velocity;          % linearized velocity based on the linear distance from each well
statematrix.referenceWell = referencewell;      % reference well for distance along each trajectory
statematrix.linearDistanceToWells = lindist;    % linearized distance to each well
statematrix.nonStandardNumSegs = nonstandard;   % the number of segments deviated off the defined trajectory; sum(statematrix.nonStandardNumSegs>0) should equal sum(statematrix.nonstandardSegmentFlag>0) calculated below


[statematrix.traj, statematrix.lindist, statematrix.nonstandardSegmentFlag, statematrix.nonstandardDirFlag] = gettraj(statematrix, segmenttable, trajwells);

statematrix.headdir = headdir;
end

% END OF MAIN FUNCTION
%---------------------------------------------------------------------

function [segmenttable, wellindex, trajwells, wellSegmentInfo, segmentInfo] = getSegmentTable(firstcoord)
% make a lookup table where each row is [traj trajsegment segID]
%
% well index gives the traj and segment index for each unique trajectory
% endpoint (well location).  If a well is used in multiple trajecories, only
% the index to the first found trajectory is used for that well.
%
% for each trajectory, a row entry in trajwell gives the start and end well
% indeces into wellindex
%
% wellSegmentInfo gives information specific to each well, such as the
% linear distance from that well to each segment.
%
% segmentInfo gives information about each segment, such as segment length

tmpcoord = [];
coordtraj = [];
trajlength = [];
wellindex = [];
wells = [];
trajwells = [];

%first, we want to make a list of all the unique wells and all the unique
%segments from the trajectories given
for i = 1:length(firstcoord)
     %add the first and last coordinates of the trajectory as the well
     %locations if these wells have not already been added to the list
          
     if ~rowfind(firstcoord{i}(1,:),wells)
        wellindex = [wellindex; [i 1]];
        wells = [wells; firstcoord{i}(1,:)];
        trajwells(i,1) = size(wellindex,1);
     else
         trajwells(i,1) = rowfind(firstcoord{i}(1,:),wells);
     end
     if ~rowfind(firstcoord{i}(size(firstcoord{i},1),:),wells)
        wellindex = [wellindex; [i size(firstcoord{i},1)]];
        wells = [wells; firstcoord{i}(size(firstcoord{i},1),:)];
        trajwells(i,2) = size(wellindex,1);
     else
         trajwells(i,2) = rowfind(firstcoord{i}(size(firstcoord{i},1),:),wells);
     end
          
     tmpcoord = [tmpcoord ; [firstcoord{i}(1:end-1,:) firstcoord{i}(2:end,:)]]; %get the coordinates of the start and end points of each segment
     tmpcoordtraj = ones(size(firstcoord{i},1)-1,1)*i;
     tmpcoordseg = (1:(size(firstcoord{i},1)-1))';
     trajlength(i,1:2) = [size(firstcoord{i},1)-1 size(coordtraj,1)];
     coordtraj = [coordtraj;[tmpcoordtraj tmpcoordseg]];
end

%any repeating segments are double labeled instead of given a unique
%segment number.
segmentCoords = [];
coordind = [];
indcount = 0;
for i = 1:size(tmpcoord,1)
    if ~isempty(segmentCoords)  
    [test, testind] = ismember(tmpcoord(i,:),segmentCoords,'rows');
    else
       test = zeros(size(tmpcoord(i,:),1),1);
       testind = 0;
    end
    %%%%%%%%%%%%%%%%%%%%%%
    
    
    reverse = 0;
    %if the segment wasn't recorded in one direction, check the reverse
    if ~(test)
        if ~isempty(segmentCoords)
            [test, testind] = ismember([tmpcoord(i,3:4) tmpcoord(i,1:2)],segmentCoords,'rows');
            reverse = 1;
        else
            test = zeros(size(tmpcoord(i,:),1),1);
            testind = 0;
        end
    end
    origin = tmpcoord(i,1:2);
    endpoint = tmpcoord(i,3:4);
    seglength = sqrt( ((endpoint(1) - origin(1))^2) + ((endpoint(2) - origin(2))^2) );
    if ~(test)
        %this segment hasn't been added to the new list yet, so give it a
        %new index
        segmentCoords = [segmentCoords; tmpcoord(i,:)];
        indcount = indcount+1;
        segmentLength(indcount) = seglength;
        coordind = [coordind; indcount];
        %record which trajectories this segment belongs to
        segmentTrajectories(indcount,1) = bitset(0,coordtraj(i,1),1);                
    else
        %this segment is already in the new list, so give it the original
        %index
        coordind = [coordind;testind];
        segmentTrajectories(testind) = bitset(segmentTrajectories(testind),coordtraj(i,1),1);
    end
end
segmentInfo.segmentCoords = segmentCoords;
segmentInfo.segmentLength = segmentLength;
segmentInfo.segmentTraj = segmentTrajectories;
segmenttable = [coordtraj coordind];

%find which segments connect to the start and end of each segment
for i = 1:size(segmentCoords,1)
    tmp = find( ((segmentCoords(:,1) == segmentCoords(i,1))&(segmentCoords(:,2) == segmentCoords(i,2))) | ((segmentCoords(:,3) == segmentCoords(i,1))&(segmentCoords(:,4) == segmentCoords(i,2))) );
    startLinkSegments{i} = setdiff(tmp,i);   
    tmp = find( ((segmentCoords(:,1) == segmentCoords(i,3))&(segmentCoords(:,2) == segmentCoords(i,4))) | ((segmentCoords(:,3) == segmentCoords(i,3))&(segmentCoords(:,4) == segmentCoords(i,4))) );
    endLinkSegments{i} = setdiff(tmp,i);  
end

%Next we want to create a connectivity matrix stating which segments are
%directly connected.  This will be used to calculate the distance from any segment to any other segment.
%Each value in the matrix is defined as e^(-distance) from segA to segB, or
%e^(-lengthSegA).  This allows us to multiply the values and get summed distance. 
connectivityTable = [];
wellsegments = [];
uniqueSegments = unique(segmenttable(:,3))';
distanceDivisor = 1000;
connectedSegments = {};
if (length(uniqueSegments) > 1)
    for i = uniqueSegments
        segindex = find(segmenttable(:,3) == i)'; % find all instances of this segment
        for j = segindex
            traj = segmenttable(j,1);
            segnum = segmenttable(j,2); % number of segment within this trajectory
            seglength = segmentLength(i);
            %find the segment indeces of the trajectory endpoints
            if (sum((wellindex(:,1) == traj)&((wellindex(:,2) == segnum)|(wellindex(:,2) == (segnum+1)))))
                transindex = min(find((wellindex(:,1) == traj)&((wellindex(:,2) == segnum)|(wellindex(:,2) == (segnum+1)))));
                wellsegments = [wellsegments; [i transindex]];
            end

            %find the segment indeces of all segments that are directly
            %connected to this segment             
            thisSegCoord = segmentCoords(i,:);
            connectedSegments_tmp = find( ((segmentCoords(:,1) == thisSegCoord(1))&(segmentCoords(:,2) == thisSegCoord(2))) | ...
                                       ((segmentCoords(:,1) == thisSegCoord(3))&(segmentCoords(:,2) == thisSegCoord(4))) | ...
                                       ((segmentCoords(:,3) == thisSegCoord(1))&(segmentCoords(:,4) == thisSegCoord(2))) | ...
                                       ((segmentCoords(:,3) == thisSegCoord(3))&(segmentCoords(:,4) == thisSegCoord(4))) ); 
            
            connectedSegments_tmp = unique(setdiff(connectedSegments_tmp,i));
            %connectedSegments = segmenttable(find((segmenttable(:,1) == traj) & ( (segmenttable(:,2) == segnum+1)|(segmenttable(:,2) == segnum-1) )),3);
            
            distanceDivisor = 1000;
            for k = connectedSegments_tmp'
                %the (i,k)th entry in the connectivity matrix says if the
                %ith segment is directly connected to the kth segment,
                %and by what distance: (exp(-lengthOfSegmentI))
                connectivityTable(i,k) = exp(-seglength/distanceDivisor); 
            end
        end
        % <traj#> cells, in which each cell contains the segments directly connected to the segment <traj#>
        connectedSegments{i} = connectedSegments_tmp;
    end
else
    %there is only one segment
    wellsegments = [1 1];
    connectivityTable = exp(-segmenttable(1,3));
end
wellsegments = sortrows(wellsegments,2);
wellsegments = wellsegments(:,1);


for i = 1:length(wellsegments)
    %calculate the length of the arms leading up to the wells (until an
    %intersection is hit). This will be used later to determine if the
    %animal completed a trajectory
    armlength(i) = 0;
    segcount = 1;
    foundIntersection = 0;
    segindex = wellsegments(i);
    tempTable = connectivityTable;
    while( (foundIntersection == 0) & (segcount < length(uniqueSegments)) )
        if(sum(tempTable(segindex,:) > 0) > 0.9)  % changed to 0.9 from 1, 1/9/17
            foundIntersection = 1;            
        else
            %we have not reached an intersection yet, so find the next
            %segments in the connection tree
            tempTable = tempTable*connectivityTable;
            segcount = segcount + 1;
        end
    end
    %calculate the length of the independent arm (all nonzero numbers on
    %the row should be the same, so just pick the maximum one)           
    armlength(i) = -log(max(tempTable(segindex,:)));      
end

%calculate the linear length from each segment to every other segment
%and the sequence of segments to get from segA to segB
tempTable = connectivityTable;
distanceTable = connectivityTable;
pathTable = cell(size(connectivityTable));
diagIndeces = find(eye(size(connectivityTable,1))); %these are the indeces to the diagonal entries
[nonZerosi,nonZerosj] = find(distanceTable > 0);
for i = 1:length(nonZerosi)
    pathTable{nonZerosi(i),nonZerosj(i)} = [nonZerosi(i) nonZerosj(i)];
end
    
while(sum(distanceTable(:) == 0) > 0)   %we do this loop until we have found a path between every pair of segments 
    
    
    oldtempTable = tempTable;
    tempTable = tempTable*connectivityTable;  %this finds the grandchildren...greatgrandchildren...and so on, of each segment, and the corresponding distance 
    zerovals = find(distanceTable == 0);
    %which segment pairs became connected?
    [switchedValsi, switchedValsj] = find( (distanceTable == 0) & (tempTable > 0) );
    
    %we only update the values that were zero (because we already found the
    %minimum path for the other pairs)
    distanceTable(zerovals) = tempTable(zerovals);
    % creates a table of every possible path between n connected segments
    for i = 1:length(switchedValsi)
        if (switchedValsi(i) ~= switchedValsj(i))
            %find which segment links one path to the new segment
            
            linker = find( (oldtempTable(switchedValsi(i),:) > 0) & (connectivityTable(:,switchedValsj(i))' > 0) );
            prepath = pathTable{switchedValsi(i),linker}(1:end-1);
            postpath = pathTable{linker,switchedValsj(i)}(2:end);
            %add the path to pathTable
            pathTable{switchedValsi(i),switchedValsj(i)} = [prepath linker postpath];
        end
    end 
end
for i = 1:length(diagIndeces)
    %the diagonal entries are wrong, so we fix them
    distanceTable(diagIndeces(i)) = 1;
    pathTable{diagIndeces(i)} = i;
end
    
%convert the exponent distance back to normal distance
%distanceTable =(-log(distanceTable)) - log(multtable);
distanceTable =-log(distanceTable) * distanceDivisor; 

%only keep the distances from the wells (we don't care about the other
%segment pairs)
distanceTable = distanceTable(wellsegments,:);
pathTable = pathTable(wellsegments,:);
segmentDirection = [];

%finally, we need to calculate whether a segment is aligned in the foreward
%or backward direction relative to each well
for i = 1:size(pathTable,1)
    for j = 1:size(pathTable,2)        
        if ( (isempty(startLinkSegments{j})) & (isempty(endLinkSegments{j})) )
            %there is only one segment
            if (i==1)
                segmentDirection(i,j) = 1;
            else
                segmentDirection(i,j) = 0;
            end              
        else
            foundpath = 0;
            for k = pathTable{i,j}
                if (ismember(k,startLinkSegments{j}))
                    segmentDirection(i,j) = 1; %from this well, the segment will linearize in the foreward direction
                    foundpath = 1;
                    break;
                elseif (ismember(k,endLinkSegments{j}))
                    segmentDirection(i,j) = 0; %or from the backward direction
                    foundpath = 1;
                    break;
                end
            end
            if (foundpath == 0)
                %this is the same segment as the well, and it will
                %therefore have either no start segments or no end
                %segments
                if (isempty(startLinkSegments{j}))
                    segmentDirection(i,j) = 1;
                elseif (isempty(endLinkSegments{j}))
                    segmentDirection(i,j) = 0;
                end
            end            
        end
    end
end

%create the wellSegmentInfo structure          
wellSegmentInfo.distanceTable = distanceTable;
wellSegmentInfo.segmentIndex = wellsegments;
wellSegmentInfo.distanceToIntersection = armlength*1000;  %convert to cm
wellSegmentInfo.wellCoord = wells;
wellSegmentInfo.pathTable = pathTable;
wellSegmentInfo.segmentDirection = segmentDirection;   

segmentInfo.connectedSegments = connectedSegments; % MS added 1/9/17 - will be used to look for deviations from defined trajectories onto neighboringh segments
end

% -------------------------------------------------------
function trajmatrix = gettraject(pos, welllocations,trajwells, lindist, wellSegmentInfo, welldist, mindiff, welldistthresh,branchpointdist);
%creates a 3 column matrix, where each row descibes a complete well-to-well
%trajectory. The columns are [starttime startwell endwell]
%welldist is the detection radius around each well (in cm)
%mindiff is the minimum time (in seconds) between detection to be called a trajectory
%mindist (welldistthresh or branchpointdist) is the threshold distance (in cm) for the program to add a trajectory when the animal turns around
%before reaching the end well

nfoodwells = size(welllocations,1);
exittimes = [];
minsamples = round(mindiff/(pos(2,1)-pos(1,1)));

% find the times that the animal was within welldist of the foodwells
for i = 1:nfoodwells
    
    tmpwellloc = squeeze(welllocations(i,:,:))';
    % find the times where the animal exited the bounding circle around the well.
    etimes = pos(find(dist(tmpwellloc, pos(:,2:3)) < welldist),1);
    
    if (~isempty(etimes))
        
        % remove all points less than mindiff apart in time
        % this will leave us with one time point per exitpoint indicating when
        % the animal left that exitpoint on each trajectory
        tmpdiff = etimes(2:length(etimes)) - etimes(1:length(etimes) - 1);
        % the leaving times for this exitpoint are the first time of each pair whose
        % difference is larger than mindiff;
        tmp = find(tmpdiff > mindiff);
        
        % tmp omits the last exitpoint, add in the last etime
        tmparray = zeros(length(tmp)+1,2);

        tmparray(:,1) = [etimes(tmp) ; etimes(end)];
        tmparray(:,2) = i;
        exittimes = [exittimes ; tmparray];
        
    end
end

if ~isempty(exittimes)
    exittimes = sortrows(exittimes, 1);
end

% this loop detects if the animal went down part more than mindist of an unnoticed arm.
% if so, then a trajectory is added. Do do this, I use a loop that breaks
% and starts over every time a new trajecory is added.  This ensures that any other arm
% entries during that period are also detected.

foundone = 0;
stopcue = 0;

while(stopcue == 0)
    foundone2 = 0;
    for ct = 2:size(exittimes,1)
        tmptraj = [exittimes(ct-1,2) exittimes(ct,2)];
        tmptimes = [exittimes(ct-1,1) exittimes(ct,1)];
        tmpindex = find((pos(:,1)> tmptimes(1)) & (pos(:,1) < tmptimes(2)));
        if (~isempty(tmpindex))
            for checkout = 1:size(lindist,2) %iterate through wells
                if ~ismember(checkout, tmptraj)
                    if isfinite(branchpointdist)
                        %use the distance from the branch point: is the linear distance from the well less than the
                        %distance from that well to its intersection, minus the branchpointdist
                        inarmind = find( (lindist(tmpindex,checkout) < (wellSegmentInfo.distanceToIntersection(checkout))-branchpointdist));
                    else
                        %use the threshhold distance from the well
                        inarmind = find( (lindist(tmpindex,checkout) < (wellSegmentInfo.distanceToIntersection(checkout))) & ...
                                     (lindist(tmpindex,checkout) < welldistthresh) );
                    end                    
                    if ~isempty(inarmind)
                        if (length(inarmind) > 1)
                            [tmpmin,endind] = min(find((diff(inarmind)>=minsamples)));
                            if isempty(endind)
                                endind = length(inarmind);
                            end
                        else
                            endind = 1;
                        end
                        [minval, minind] = min(lindist(tmpindex(inarmind(1:endind)),checkout));
                        outtime = pos(tmpindex(inarmind(minind)),1);
                        exittimes = [exittimes(1:ct-1,:);outtime checkout;exittimes(ct:end,:)];
                        foundone = 1;
                        foundone2 = 1;                       
                        break;
                    end
                end
            end          
        end
        if (foundone)
            %a trajectory was added, so we break and start over
            foundone = 0;
            break;
        end
    end
    if (foundone2 == 0)
        %no trajectory was added, so we are done
        stopcue = 1;
    else
        foundone2 = 0;
    end
end

%create a matrix describing the start time, start well, and end well of
%each complete trajectory
trajmatrix = [];
for i = 1:size(exittimes,1)-1
    trajmatrix(i,1:3) = [exittimes(i,1) exittimes(i,2) exittimes(i+1,2)];
end


end


%------------------------------------------------------------------
function [state, lindist, nonstandard_segment, nonstandard_headdir] = gettraj(statematrix, segmenttable, trajwells)

%get the most probable linear trajectory number, based on well exit/enters
%and head direction.

traject = trajwells;
includeStates = 1:6;

%if no reference well is assigned, it is not a valid point
%also, only points inside the designated time range are valid
validPoints = ((statematrix.referenceWell > 0)); 

validPointsIndex = find(validPoints); 
trajvector = ones(length(validPointsIndex),1)*-1;
statevector = ones(length(validPointsIndex),1)*-1;
referencewell = statematrix.referenceWell(validPointsIndex);
welltraj = statematrix.wellExitEnter(validPointsIndex,:); %exit and enter wells
segmentindex = statematrix.segmentIndex(validPointsIndex);
samewell = (welltraj(:,1) == welltraj(:,2)); %Which times are during trajectories between the same wells 
diffwell = (welltraj(:,1) ~= welltraj(:,2)); %Which times are during trajectories between different wells
trajcount = rowcount(statematrix.segmentIndex(validPointsIndex),segmenttable(:,3)); %how many of the linear trajectories does each segment belong to

matrixIndex = sub2ind(size(statematrix.linearDistanceToWells),validPointsIndex,referencewell);  % simply the indices for valid ponits to the REFERENCE (well 1) well
lineardistance = statematrix.linearDistanceToWells(matrixIndex);  % picks out the linear distance to the REFERENCE WELL (given matrixIndex)

forwarddir = ((statematrix.segmentHeadDirection(matrixIndex) >= 0));       %facing positve dir and moving positive dir (both head and motion
                                                                            %direction in the positive direction, and velocity greater than minvelocity)
backwarddir = ((statematrix.segmentHeadDirection(matrixIndex) < 0));        %facing negative dir and moving negative dir

alltraj = [];
%compile the full list of trajectories
for trajnum = 1:size(traject,1)
    currtraj = traject(trajnum,:);
    alltraj = [alltraj; currtraj 2*(trajnum)-1];
    alltraj = [alltraj; currtraj([2 1]) 2*(trajnum)];
end

%convert segmenttable to a version labeled by alltraj (i.e. 1:30 for multiple W)
new_segmenttable = sortrows([horzcat(segmenttable(:,1)*2-1,segmenttable(:,2:3));horzcat(segmenttable(:,1)*2,segmenttable(:,2:3))],1); 

%LEVEL 1   
for trajnum = 1:size(alltraj,1)
    currtraj = alltraj(trajnum,1:2);
    
    %which times occur when the animal is intransit between the two wells
    %for the current trajectory
                                   % exit well == start well  enter well == end well of this traj   %can now exclude because already separated by direction:  exit well == end well of this traj   enter well == start well                      
    inDefinedTraj(:,trajnum) = ((welltraj(:,1)==currtraj(1))&(welltraj(:,2)==currtraj(2))); % | ((welltraj(:,1)==currtraj(2))&(welltraj(:,2)==currtraj(1))));   
        % inDefinedTraj
            % 1st col: traj 1  // 2nd col: traj 2
            % row corresponds to each time point
        
    %which segments are in the current trajectory
    segmentsInTraj = new_segmenttable(find(new_segmenttable(:,1) == trajnum),3);
    
    %find all times when the animal was on the current trajectory (either
    %in the foreward or backward direction).  This is defined as the times 
    %when the animal is in transit between the two correct endpoints and is
    %on a track segment that is part of the trajectory.
    
    %first select the direction 
    if floor(trajnum/2)==trajnum/2 %is even, going in backward direction
        usedir = backwarddir;
    elseif floor(trajnum/2)~=trajnum/2 %is odd, going in forward direction
        usedir = forwarddir;
    end
    
     % Here searching for statevec indices such that  <times in current trajectory> & <times in current trajectory's segments> & <movement in the chosen direction>
    findindex = find( inDefinedTraj(:,trajnum) & ismember(segmentindex,segmentsInTraj) & usedir);   
    %assign the current trajectory to those indeces
    %Note: if direction does not match, the index will currently be left unassigned
    trajvector(findindex) = trajnum;         % odd outputs -- moving from left well to right well (in pos recon view)
                                             % even outputs -- moving from right well to left well

    if ismember(1,includeStates)
        statevector(findindex) = trajvector(findindex);
    end
    
end

undefinedindex = (trajvector == -1); %still undefined times
inAnyDefinedTraj = (sum(inDefinedTraj,2) > 0); %These times occur when the animal is in transit between two wells of a defined trajectory


%LEVEL 2
%for the times that are still undefined :
%   1) not occurring during any of the above well-to-well trajectories, 
%   2) occur when the animal is not leaving and entering the same well,
%   3) is in a track segment that only belongs to one of the defined linear trajectories, assign the proper trajectory 
findindex = find((undefinedindex) & (diffwell) & (trajcount == 1) & (forwarddir) & (~inAnyDefinedTraj)); 
% in multi-W, no segment belongs to only one traj, so this will always be empty
if ~isempty(findindex)
    trajvector(findindex) = 2*(segmenttable(rowfind(segmentindex(findindex),segmenttable(:,3)),1))-1;
    if ismember(2,includeStates)
        statevector(findindex) = trajvector(findindex);
    end
    findindex = find((undefinedindex) & (diffwell) & (trajcount == 1) & (backwarddir) & (~inAnyDefinedTraj));
    trajvector(findindex) = 2*(segmenttable(rowfind(segmentindex(findindex),segmenttable(:,3)),1));
    if ismember(2,includeStates)
        statevector(findindex) = trajvector(findindex);
    end
    
    if ((nargout == 1) & (max(includeStates) <= 2))
        return
    end
end



%LEVEL 3
%fill all times when the animal is on a segment that belongs to only one
%trajectory (even if it is leaving and entering the same well)
% again for multi-W, this will always be empty
% FORWARD
undefinedindex = (trajvector == -1); %still undefined times
findindex = find((undefinedindex) & (trajcount == 1) & (forwarddir)); 
if ~isempty(findindex)
    % this next line: 1) looks up which segment (in particular, segment 1-5 index) the animal is currently on
    % 2) from this segment, looks up the linear trajectory (1,2)
    % 1: to the RIGHT  >> assigns trajvector = 1
    % 2: to the LEFT   >> assigns trajvector = 3
    trajvector(findindex) = 2*(segmenttable(rowfind(segmentindex(findindex),segmenttable(:,3)),1))-1;
    if ismember(3,includeStates)
        statevector(findindex) = trajvector(findindex);
    end
    % BACKWARD
    findindex = find((undefinedindex) & (trajcount == 1) & (backwarddir));
    % this next line: 1) looks up which segment (in particular, segment 1-5 index) the animal is currently on
    % 2) from this segment, looks up the linear trajectory (1,2)
    % 1: to the RIGHT  >> assigns trajvector = 2
    % 2: to the LEFT   >> assigns trajvector = 4
    trajvector(findindex) = 2*(segmenttable(rowfind(segmentindex(findindex),segmenttable(:,3)),1));
    if ismember(3,includeStates)
        statevector(findindex) = trajvector(findindex);
    end
end

% LEVEL 4 (MS added: important for mulit-W)
% Fill in the undefined times when the animal was on an ambiguous segment (like the
%home arm in a w-maze) which belongs to multiple traj, and startwell is the
%same as the endwell (trackback) --> these get assigned their own
%trajectory numbers
undefinedindex = (trajvector == -1); %the undefined times 
findindex = find((undefinedindex) & (samewell) & (trajcount > 1)); 
trackbackIDs = 31:36;
%trackbacks = welltraj(findindex,:);
% iterate through wells
for w = unique(trajwells)'
    iswell = (welltraj(:,1)==w & welltraj(:,2)==w); % find times that start and end with this well
    findex2 = findindex(ismember(findindex,find(iswell))); % find indices that match
    trajvector(findex2) = trackbackIDs(w);
    statevector(findex2) = trajvector(findex2);
end

% at this point, undefined indices might still result from forwarddir or
% backwarddir == 0, which means headdir was probably NaN because all points are already valid 

%LEVEL 5 (MS added: important for mulit-W)
%Fill in the undefined times when the animal was on an ambiguous segment and there is no defined direction OR there was a 
% head direction that doesn't match the trajectory (e.g. brief diode occlusion made the headdir look negative on a forward traj)
nonst_segment = zeros(size(trajvector)); %binary vector to indicate whether a traj was smooth/standard (0) or nonstandard (1)
nonst_headdir = zeros(size(trajvector)); %binary vector to indicate whether a headdirection was nonstandard (1) compared to the trajectory (forward or backward dir)
possibletraj = zeros(length(trajvector),2);
trajToMerge = zeros(size(trajvector));
undefinedindex = (trajvector == -1); %the undefined times 

possibletraj(undefinedindex,:) = welltraj(undefinedindex,:);
trajToMerge(undefinedindex) = rowfind(possibletraj(undefinedindex,:),alltraj(:,1:2));

% iterate through possible traj, ask if the segments at these indices are
% part of segmentsInTraj: if not, merge with the current traj as indicated by well start and end but flag as "nonstandard" 
for trajnum = unique(trajToMerge(trajToMerge~=0))'
    segmentsInTraj = new_segmenttable(find(new_segmenttable(:,1) == trajnum),3);
    findindex = find((undefinedindex) & (trajToMerge==trajnum) & ~ismember(segmentindex,segmentsInTraj));
    findindex2 = find((undefinedindex) & (trajToMerge==trajnum) & ismember(segmentindex,segmentsInTraj));
    trajvector(findindex) = trajToMerge(findindex); %because of atypical segment visited
    trajvector(findindex2) = trajToMerge(findindex2); %often because of segmentHeadDirection NaNs or spurious changes in head direction against the direction of the trajectory
    nonst_segment(findindex) = 1;
    nonst_headdir(findindex2) = 1;
    statevector(findindex) = trajvector(findindex);
    statevector(findindex2) = trajvector(findindex2);
end


% Check for brief changes in head direction that got erroneously read out
% as a change in trajectory

% trajchanges = [NaN; diff(statevector)]; %pad w/ NaN so it will be the same length
% changeinds = find(trajchanges~=0);
% trajlengths = [NaN; diff(changeinds)];


%%%%% DEPRECATED BELOW HERE: MS left commented for reference
%if facing positive direction, then assign the traj of the next traj in the
%future
% findindex = find((undefinedindex) & (forewarddir)); 
% trajvector = vectorfill(trajvector, -1, 1, findindex);
% findex2 = findindex(find(~mod(trajvector(findindex),2)));
% trajvector(findex2) = trajvector(findex2)-1;
% if ismember(4,includeStates)
%     statevector(findindex) = trajvector(findindex);
%     nonstandard(findindex) = 1; % set to 1 to flag as a "nonstandard" trajectory, there was a diversion off the direct course
% end
% %otherwise assign the closest traj that happened in the past 
% findindex = find((undefinedindex) & (trajcount > 1) & (backwarddir)); 
% trajvector = vectorfill(trajvector, -1, -1, findindex);
% findex2 = findindex(find(mod(trajvector(findindex),2)));
% trajvector(findex2) = trajvector(findex2)+1;
% if ismember(4,includeStates)
%     statevector(findindex) = trajvector(findindex);
% end
% 
% 
% 
% %LEVEL 6
% % now we may have some undefined points in the beginning and end of the
% % session- we assign these to the first trajectory
% undefinedindex = (trajvector < 0); %still undefined times
% trajvector(find(undefinedindex)) = -1;
% statevector(find(undefinedindex)) = -1;
% findindex = find((undefinedindex) & (forewarddir)); 
% trajvector(findindex) = 1;
% if ismember(5,includeStates)
%     statevector(findindex) = trajvector(findindex);
% end
% findindex = find((undefinedindex) & (backwarddir)); 
% trajvector(findindex) = 2;
% if ismember(5,includeStates)
%     statevector(findindex) = trajvector(findindex);
% end
% 
% 
% 
% %LEVEL 7
% % Fill in all remaining points, regardless of running velocity
% undefinedindex = find(trajvector == -1); %still undefined times
% trajvector = vectorfill(trajvector, -1, 0, undefinedindex);
% if ismember(6,includeStates)
%     statevector = trajvector;
% end




state = ones(length(statematrix.time),1)*-1;
lindist = ones(length(statematrix.time),1)*-1;
nonstandard_segment = ones(length(statematrix.time),1)*-1;
nonstandard_headdir = ones(length(statematrix.time),1)*-1;
lindist = ones(length(statematrix.time),1)*-1;
state(validPointsIndex) = statevector;
lindist(validPointsIndex) = lineardistance;
nonstandard_segment(validPointsIndex) = nonst_segment;
nonstandard_headdir(validPointsIndex) = nonst_headdir;
end



