function coords = getcoord_6armtrack(directoryname,pos,index,vidframe,cmperpix)
%edited 11-14-06 by AS
%updated 01-05-17 by MS
%updated 01-15-17 by DR
%this program is called by CREATETASKSTRUCT to produce the trajectory
%coodinates for a wtrack.
%Click locations in the following order:
%   
%   1    2    3    4    5    6
%   |    |    |    |    |    |
%   |    |    |    |    |    |
%   |    |    |    |    |    |
%   |    |    |    |    |    |
%   7----8----9---10---11----12
%



fid = figure;
if exist('vidframe','var')
    %plot image and flip l/r
    image(vidframe)
    set(gca,'YDir','normal') %flips both image and axis u/d to align with pos recon coordinates
    hold on
    plot((pos(:,2)/cmperpix(index(2))),(pos(:,3)/cmperpix(index(2))), ':'); 
else
    plot(pos(:,2),pos(:,3));
end

[x,y] = ginput(12);
%%%% Must order the first well of each trajectory so they get assigned 1-6
% from left to right (also remember axis is flipped vertically)
%%%% This is all possible trajectories in one direction
%%%% Convention is to order the trajectories from lowest # well to highest
lincoord{1} = [x([1 7 8 2]) y([1 7 8 2])];      % 1 to 2 %no skips
lincoord{2} = [x([2 8 9 3]) y([2 8 9 3])];      % 2 to 3 
lincoord{3} = [x([3 9 10 4]) y([3 9 10 4])];    % 3 to 4
lincoord{4} = [x([4 10 11 5]) y([4 10 11 5])];  % 4 to 5
lincoord{5} = [x([5 11 12 6]) y([5 11 12 6])];  % 5 to 6
lincoord{6} = [x([1 7 8 9 3]) y([1 7 8 9 3])];  % 1 to 3 %skip 1
lincoord{7} = [x([2 8 9 10 4]) y([2 8 9 10 4])];% 2 to 4
lincoord{8} = [x([3 9 10 11 5]) y([3 9 10 11 5])]; % 3 to 5
lincoord{9} = [x([4 10 11 12 6]) y([4 10 11 12 6])]; % 4 to 6
lincoord{10} = [x([3 9 10 11 12 6]) y([3 9 10 11 12 6])]; % 3 to 6 (home to outer) % skip 2
lincoord{11} = [x([1 7 8 9 10 4]) y([1 7 8 9 10 4])]; % 1 to 4 (outer to home)
lincoord{12} = [x([2 8 9 10 11 5]) y([2 8 9 10 11 5])]; % 2 to 5
lincoord{13} = [x([1 7 8 9 10 11 5]) y([1 7 8 9 10 11 5])]; % 1 to 5  % skip 3
lincoord{14} = [x([2 8 9 10 11 12 6]) y([2 8 9 10 11 12 6])]; % 2 to 6  % skip 3
lincoord{15} = [x([1 7 8 9 10 11 12 6]) y([1 7 8 9 10 11 12 6])]; % 1 to 6 (outer to outer, rat actually does this a lot)

numtimes = size(pos,1);
for i = 1:length(lincoord)
    if exist('M','var')
    coords{i} = repmat(lincoord{i}*cmperpix,[1 1 numtimes]);
    else
        coords{i} = repmat(lincoord{i},[1 1 numtimes]);
    end
end
close(fid);