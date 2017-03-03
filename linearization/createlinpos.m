% [linpos, lincoordlist] = LINEARIZEPOS(posstruct, taskstruct, index)
%   01.18.17: renamed to createlinpos for clarity
%          Takes each position from the specified element of posstruct and 
%  	   projects it onto the linear trajectory segments for the environment
%  	   and task defined in the specified element of taskstruct.
%          This routine finds the segment closest to each position, projects
%          the position to the closest location on the segment. 
%	   linpos is the list of linearized positions and 
%          lincoordlist is the number of the coordinate list and the segment
%          the animal is on at each time point.

% MS modified 1.10.17 from Annabelle Singer original

function [linpos, segmentIndex, segdist] = createlinpos(pos, task, index, coord, varargin)

disp('Getting linpos...')
% pos = posstruct{index(1)}{index(2)}.data;
% task = taskstruct{index(1)}{index(2)};
linpos = pos;

z = find(pos(:,2) == 0);
nonz = find(pos(:,2));
pos = pos(nonz,:);

coordlist = zeros(size(pos,1), 2);

newpos = pos;

maxv = 500;

%set variable options
for option = 1:2:length(varargin)-1
    switch varargin{option}
        case 'maxv'
            maxv = varargin{option+1};
    end
end


% set the variable for the list of projection coordinates
trajcoord = task.linearcoord;

% project each position point onto each segment
ntraj = length(trajcoord);
nseg = length(coord);
maxsegdiff = cell(nseg,1); 

% maxsegdiff is the maximum allowed number of segments that
% could be traversed from one time point to the next
for j = 1:length(coord)
    if (size(coord{j},1) <= 3)
	% there are only one or two segments, so maxsegdiff isn't really 
	% meaningful. We set it to 10 so that it is ignored
	maxsegdiff{j} = 10;
    else
	maxsegdiff{j} = min([(length(coord{j})-2) maxsegdiff{j}]);
    end
end
inbound = 0;
lastvalid = -1;
segdist = zeros(size(pos,1),1);
for i = 1:size(pos,1)
    % project each position point onto the linear segments
    if (pos(i,2) == 0)
        % this is an invalid point, so set newpos to zero for this position
        newpos(i, 2:3) = 0;
    else
        tmppos = [];
	for j = 1:nseg
        
        % here input to projectpoint a SINGLE (x,y) coord (current 2D position), and the two (x,y) endpoint coords of a line segment
            try    % using smoothed position, if it exists      
            tmppos{j} = projectpoint(pos(i,6:7), coord{j});
            % tmppos :
                % rows: one row for each segment
                % [sm-x sm-y distance onseg segnum] =
                % [xpos-of-closest-point-on-seg, ypos-of-closest-point-on-seg, distance-to-closest-point, whether-projected-point-was-within-segment, segment-number]
            catch
%                 if do_once == 0
                disp('(!!) using x and y from pos{}{} instead of sm-x and sm-y..')
%                 do_once = 1;
%                 end
                tmppos{j} = projectpoint(pos(i,2:3), coord{j});
            end

        % change the segment number to equal the segment index
        tmppos{j}(:,5) = j;
	end
        % take the point that is the least distance from the segments and that
        % does not force an excessive velocity from the last point
        % also check to see if the closest point forces a transition between
        % the first and the last segment
        newind = cell(length(coord),1);
        if (lastvalid ~= -1)
            for j = 1:length(tmppos)
                if (~isempty(tmppos{j}))
                    tmppos{j} = sortrows(tmppos{j},3);
                    % check the velocity to the last valid point
                    tmppnt = zeros(size(tmppos{j},1),2);
                    tmppnt(:,1) = newpos(lastvalid,2);
                    tmppnt(:,2) = newpos(lastvalid,3);
                    vel = dist(tmppos{j}(:,1:2), tmppnt) ./ (newpos(i,1) - ...
                        newpos(lastvalid,1));
                    % the next point is the first point in the list with a
                    % corresponding velocity less than maxv
                    segdiff = Inf;
%		if pos(i,2) > 200
%		    keyboard
%		end
                    newind{j} = min(find(vel < maxv));
                    %newind{j} = 1;
                    currentseg = tmppos{j}(newind{j}, 5);
                    % if the closest point forces a transition between the
                    % first and the last segment, use the next closest point.
                    % This works both within and across linearcoord lists
                    if ((abs(currentseg - lastseg)) == maxsegdiff{j})
			sprintf('correcting segment on point %d', i)
                        ind = find(vel < maxv);
                        if (size(ind,1) > 1)
                            newind{j} = ind(2);
                        else
                            newind{j} = [];
                        end
                    end
                end
            end
            % newind{j} contains the closest point on the jth linearcoord
            % list, so compare the distances of the points from each
            % linearcoord list to choose the closest point.
	    empty = 1;
	    for j = 1:nseg %ntraj
		if (~isempty(newind{j}))
		    empty = 0;
		end
	    end
	    if (empty)
                % set this point to zero, as it violates the velocity cutoff for
                % both linearcoord lists or the original point was zero
		sprintf('Warning: zero position at position element %d, terpos{%d}{%d}\n', i, index(1), index(2))
		
                newpos(i,2:3) = [0 0];
            else
                clist = 0;
                mindist = Inf;
                for j = 1:length(newind)
                    if ((~isempty(newind{j})) & (tmppos{j}(newind{j},3) < ...
                                                  mindist))
                        mindist = tmppos{j}(newind{j},3);
                        clist = j;
			% save the coordinate list number and the segment
			% number
			coordlist(i,:) = [clist tmppos{j}(newind{j},5)];
                    end
                end
                currentclist = clist;
                currentseg = tmppos{clist}(newind{clist},5);
                newpos(i,2:3) = round(tmppos{clist}(newind{clist}, 1:2));
                lastseg = tmppos{clist}(newind{clist}, 5);
                segdist(i) = sqrt(sum((tmppos{clist}(1:2) - coord{clist}(1,:)).^2)); %the distance along the segment  
                lastvalid = i;
            end
        else 
            % the last point was not valid so pick the point with the minimum 
            % distance
            tmp = [];
            for j = 1:length(tmppos)
                tmp = [tmp tmppos{j}(1,3)]; 
            end
            [v, newind] = min(tmp);
            newpos(i,2:3) = round(tmppos{newind}(1, 1:2));
            lastseg = tmppos{newind}(1, 5);
            lastvalid = i;
        end
    end
end

% put the zeros back in, if they exist

linpos(nonz,2:3) = newpos(:,2:3);
lincoordlist = zeros(length(linpos(:,1)),2);
lincoordlist(nonz,:) = coordlist;
segmentIndex = coordlist(:,1);



