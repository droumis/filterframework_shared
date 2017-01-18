

%% Type in animal, day, and indices to probe here manually
% By Kenny Kay
% Mari updated to run through multiple W segments and trajectories
% use this script to check segment and trajectory assignment by linearizeposition

% can also be used to figure out each:
% [day epoch startposind endposind segmentnum]
% i.e, the linpos indices that need to be manually reassigned to the
% correct trajectory
% startposind -- corresponds to the pos{}{} rows

% these are then manually entered into manual_segment option of
% kk_lineardayprocess -- downstream, kk_linearizeposition recognizes these
% places to force the correct segment

% how to use: use plot_lindistseg and plot_xysegment together
% the faulty areas will be extremely clear on plot_xysegment
% find where these segments occur in the plot_lindistseg, then zoom in
% manually to find the FIRST index (and a stretch afterwards) to change
%
% set segment and trajectory colors below

animal = 'Geronimo';
%animalinfo = animaldef(animal);
animaldir = '/opt/data40/mari/Ger/'; %animalinfo{2};
animalprefix = 'ger'  %animalinfo{3};
day = 5;
indexvec_flag = 1;        % leave this toggled -- directly identifies indices
startval = 1; %ununsed (?)
endval = 2500; %unused (?)
linear_flag = 0;        % set to 1 if linear (not W) track -- then will plot trajectories instead of segments

epoch_toplot = 3;

% plot formats
plot_lindistseg = 1;

% plot positions of identified segments
plot_xysegment = 1;
% if plot_xysegment
%     epoch_toplot = 5; %find(~cellfun('isempty',linpos{day}),1);
% end

% plot identified trajectories
plot_xytraj = 1;

% plot head dir on position
plot_xyheaddir = 0;
if plot_xyheaddir
    epoch_toplot = 5;
end


%% Load pos and linpos.
pos = loaddatastruct(animaldir, animalprefix, 'pos', day);
linpos = loaddatastruct(animaldir, animalprefix, 'linpos', day);


%% Module 1:

if plot_lindistseg
    
    counter = 1;
    
    epochs = find(~cellfun('isempty',linpos{day}));
    
    H = figure;
    
    for ep = epochs
        
        time = linpos{day}{ep}.statematrix.time;
        if ~linear_flag
            segmentnum = linpos{day}{ep}.statematrix.segmentIndex;
        else
            segmentnum = linpos{day}{ep}.statematrix.traj;
        end
        % latest version of linpos should have .lindist
        if isfield(linpos{day}{ep}.statematrix,'lindist')
            lindist = linpos{day}{ep}.statematrix.lindist;
        else  % Frank, older lin pos versions
            if unique(linpos{day}{ep}.statematrix.referenceWell) ~= 1
                error('something interesting with linpos..')
            else
                lindist = linpos{day}{ep}.statematrix.linearDistanceToWells(:,1);
            end
        end
        subplot(length(epochs),1,counter)
        
        %  x vector
        if indexvec_flag
            xvec = (1:length(time))';
        else
            xvec = time;
        end
        
        % plot grey line for all positions
        plot(xvec,lindist,'linewidth',1.5,'color',[.8 .8 .8]);
        hold on
        
        % plot lindist  +  spikes
        for segno = 1:11
            inds = [];
            inds = double(segmentnum == segno);
            inds(find(inds==0))=nan;
            if segno == 1
                clr = 'k';
            elseif segno == 2
                clr = 'r';
            elseif segno == 3
                clr = 'g';
            elseif segno == 4
                clr = 'b';
            elseif segno == 5
                clr = 'm';
            elseif segno == 6
                clr = [1 0.8 0]; %yellow
            elseif segno == 7
                clr = [0.3 0.6 0]; %dark green
            elseif segno == 8
                clr = [0 0.8 1]; %cyan
            elseif segno == 9
                clr = [1 0.6 0]; %orange
            elseif segno == 10
                clr = [0.5 0 0.5]; %purple
            elseif segno == 11
                clr = [0 0.5 0.5]; %teal
            end
            % plot animal's travel lindist
            plot(xvec.*inds,lindist.*inds,'-','linewidth',3,'color',clr);
        end
        axis tight
        
        % plot starttime to endtime trace
        if 0
            if ~indexvec_flag
                startval = lookup(startval,xvec);
                endval = lookup(endval,xvec);
            end
            plot(xvec(startval:endval),lindist(startval:endval),'-','color',[.4 .4 .4]);
        end
        
        % plot title
        title([animal(1:3) ' lindist - traj, ' ' day ' num2str(day) 'epoch ' num2str(ep)],'fontsize',16,'fontweight','bold');
        counter = counter + 1;
    end
    
end

%% Module 2: Plot xy positions colored by SEGMENT NUMBER

if plot_xysegment
    
    
    
    posdata = pos{day}{epoch_toplot}.data;
    if ~linear_flag
        segmentnum = linpos{day}{epoch_toplot}.statematrix.segmentIndex;
    else
        segmentnum = linpos{day}{epoch_toplot}.statematrix.traj; 
    end
    %     % for OLDER pos -- do not have sm-x and sm-y etc.
    %     % quick fix: recopy them and send out a warning
    %     disp('this animal does NOT have sm- pos.. simply using the older x and y')
    %     posdata(:,6) = posdata(:,2);
    %     posdata(:,7) = posdata(:,3);
    
    
    
    % for each segment numbers, plot in different colors
    for segno = 1:11
        H = figure;
        % plot grey line for all positions
        plot(posdata(:,6),posdata(:,7),'linewidth',1.5,'color',[.8 .8 .8]);
        hold on
        if segno == 1
            clr = 'k';
        elseif segno == 2
            clr = 'r';
        elseif segno == 3
            clr = 'g';
        elseif segno == 4
            clr = 'b';
        elseif segno == 5
            clr = 'm';
        elseif segno == 6
            clr = [1 0.8 0]; %yellow
        elseif segno == 7
            clr = [0.3 0.6 0]; %dark green
        elseif segno == 8
            clr = [0 0.8 1]; %cyan
        elseif segno == 9
            clr = [1 0.6 0]; %orange
        elseif segno == 10
            clr = [0.5 0 0.5]; %purple
        elseif segno == 11
            clr = [0 0.5 0.5]; %teal
        end
        inds = (segmentnum == segno);
        if ~linear_flag
            scatter(posdata(inds,6),posdata(inds,7),4,clr); %,'linewidth',2,'color',);
        else
            scatter(posdata(inds,2),posdata(inds,3),'.','linewidth',2,'color',clr);
        end
     
            title([animal(1:3) ' day ' num2str(day) ' epoch ' num2str(epoch_toplot) ' segment ' num2str(segno) ' : ' num2str([startval endval])],'fontsize',16,'fontweight','bold');

    end
    
    
    % plot starttime to endtime trajectory
%     if ~indexvec_flag
%         startval = lookup(startval,posdata(:,1));
%         endval = lookup(endval,posdata(:,1));
%     end
%     plot(posdata(startval:endval,6),posdata(startval:endval,7),'.','linewidth',6,'color','y');
%     scatter(posdata(startval,6),posdata(startval,7),1000,'k','.');
    
 pause    
end


%% Module 3: Plot xy positions w/ head direction, highlighting travel between startval and endval

if plot_xyheaddir
    
    H = figure;
    
    posdata = pos{day}{epoch_toplot}.data;
    trajnum = linpos{day}{epoch_toplot}.statematrix.traj;
    
    % plot grey line for all positions
    plot(posdata(:,6),posdata(:,7),'linewidth',1.5,'color',[.8 .8 .8]);
    hold on
    
    % plot starttime to endtime trajectory
    if ~indexvec_flag
        startval = lookup(startval,posdata(:,1));
        endval = lookup(endval,posdata(:,1));
    end
    plot(posdata(startval:endval,6),posdata(startval:endval,7),'linewidth',5,'color','r');
    % plot starttime to endtime head direction following the lines
    for ind = (find(trajnum == -1))'; % plot head dir at trajectory-undefined points   %startval:endval
        xdir = 5*cos(posdata(ind,8));
        ydir = 5*sin(posdata(ind,8));
        % plot "arrow line"
        plot([posdata(ind,6)  posdata(ind,6)+xdir],...
            [posdata(ind,7)  posdata(ind,7)+ydir],'linewidth',2,'color',[.5 .5 .5]);
        % plot point at base of arrow
        scatter(posdata(ind,6),posdata(ind,7),20,'b');
        % plot point at tip of arrow
        plot(posdata(ind,6)+xdir,posdata(ind,7)+ydir,'.','markersize',25,'color',[.2 .2 .7]);
        % if the head direction data isn't there (nan) then plot a circle
        % at the position point
        if isnan(xdir) || isnan(ydir)
            scatter([posdata(ind,6)],[posdata(ind,7)],80,'k','.');
        end
    end
    scatter(posdata(startval,6),posdata(startval,7),3000,'k','.');
    title([animal(1:3) ' ' num2str(epoch_toplot) ' : ' num2str([startval endval])],'fontsize',16,'fontweight','bold');
    
end



%% Module 4: Plot xy positions colored by trajectory NUMBER

if plot_xytraj
    
%     H = figure;
    
    posdata = pos{day}{epoch_toplot}.data;

        trajnum = linpos{day}{epoch_toplot}.statematrix.traj;

    %     % for OLDER pos -- do not have sm-x and sm-y etc.
    %     % quick fix: recopy them and send out a warning
    %     disp('this animal does NOT have sm- pos.. simply using the older x and y')
    %     posdata(:,6) = posdata(:,2);
    %     posdata(:,7) = posdata(:,3);
    

    
    % -- for each trajectory numbers, plot in different colors - two
    % directions between the same pair of wells are plotted in the same color
    % -- to check for unassigned indices, set trajno to -1
    for trajno = [1:36]
%         clr = 'k'; % to check unassigned indices
        if trajno == 1 || trajno == 2
            clr = 'k';
        elseif trajno == 3 || trajno == 4
            clr = 'r';
        elseif trajno == 5 || trajno == 6
            clr = 'g';
        elseif trajno == 7 || trajno == 8
            clr = 'b';
        elseif trajno == 9 || trajno == 10
            clr = 'm';
        elseif trajno == 11 || trajno == 12
            clr = [1 0.8 0]; %yellow
        elseif trajno == 13 || trajno == 14
            clr = [0.3 0.6 0]; %dark green
        elseif trajno == 15 || trajno == 16
            clr = [0 0.8 1]; %cyan
        elseif trajno == 17 || trajno == 18
            clr = [1 0.6 0]; %orange
        elseif trajno == 19 || trajno == 20
            clr = [0.5 0 0.5]; %purple
        elseif trajno == 21 || trajno == 22
            clr = [0 0.5 0.5]; %teal
        elseif trajno == 23 || trajno == 24
            clr = [0.8 0.6 0.8]; %mauve
        elseif trajno == 25 || trajno == 26
            clr = [0.9 0.65 0.9]; %mauve
            elseif trajno == 27 || trajno == 28
            clr = [0.6 0.6 0]; %olive
            elseif trajno == 29 || trajno == 30
            clr = [0.3 0.3 0.1]; %mauve
        elseif trajno == 31 || trajno == 32 || trajno == 33 || trajno == 34 || trajno == 35 || trajno == 36 %all trackbacks
            clr = [0.3 0.3 0.3]; %dark grey
        end
        inds = (trajnum == trajno);
        % to plot time points per trajectory that occurred on segments not
        % technically included in the trajectory
%         inds = ((trajnum == trajno) & (linpos{day}{epoch_toplot}.statematrix.nonstandardSegmentFlag>0))
%         
   % plot grey line for all positions
    figure
    plot(posdata(:,6),posdata(:,7),'linewidth',1.5,'Color',[.8 .8 .8]);
    hold on
    scatter(posdata(inds,6),posdata(inds,7),4,clr);
            %plot(posdata(inds,6),posdata(inds,7),'-','linewidth',2,'color',clr);
       title([animal(1:3) ' day ' num2str(day) ' epoch ' num2str(epoch_toplot) ' traj ' num2str(trajno)],'fontsize',16,'fontweight','bold');
        
    end
    
    
%     % plot starttime to endtime trajectory
%     if ~indexvec_flag
%         startval = lookup(startval,posdata(:,1));
%         endval = lookup(endval,posdata(:,1));
%     end
%     plot(posdata(startval:endval,6),posdata(startval:endval,7),'.','linewidth',6,'color','y');
%     scatter(posdata(startval,6),posdata(startval,7),1000,'k','.');
    
    %title([animal(1:3) ' day ' num2str(day) ' epoch ' num2str(epoch_toplot) ' : ' num2str([startval endval])],'fontsize',16,'fontweight','bold');
    
end


