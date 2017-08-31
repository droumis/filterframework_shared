function [out] = dfa_getrewardtrigspiking2(index, excludeperiods, spikes, rewardinfo, eventcons, pos, varargin)
% MS 2016: adapted from dfakk_getrewardtrigspiking(index, excludeperiods, spikes, rewardinfo, varargin)
%
%   This function transcribes reward triggered spiking for every reward trigger, as well as reward triggered ripples.

% UPDATED 11.3.16 - Trials are excluded in which rat does not stay
% below maxvelocity for the duration of the window.  This matches velocity
% for reward vs. error trials. (does velocity by default)
%
%   use singlecellanal
%
%   index [day epoch tetrode cell]
%
%   out = out.index                 [D E T C], gives the identity of the cells
%         out.times                 vector of bin times for the histograms
%         out.trigmatrix
%         out.trigmatrixdescript   'well #  //  error (0) or reward (1) // inbound (10) or outbound (11)';
%         out.spikeraster           a histogram count of spikes around the trigger in 1 ms bins
%         out.spikehist             a histogram count of spikes around the trigger in larger bins
%         out.ripples               a histogram count of ripples around the trigger
%         out.spiketimes            transcribed spike times in the window around the trigger
%         out.rewardinputtimes      Nosepoke times on reward or error trials
%         out.rewardoutputtimes     Reward delivery times on reward trials,
%                                   pseudo delivery times on error trials if pseudo is specified
%         out.rewardoutcomes        Whether trials were correct or not
%         out.velocity              animal's velocity in the window around the trigger time

% default options
pseudodelay = 2;
minthresh = 3;
maxvelocity = 4;
window = [5 5];  % in sec
binsize = 0.001;  % in sec
psthbinsize = 0.01; %in sec, bins for psth
inputtimes = 1;
outputtimes = 0;
doripples = 1;

for option = 1:2:length(varargin)-1
    switch varargin{option}
        case 'window'
            window = varargin{option+1};
        case 'binsize'
            binsize = varargin{option+1};
        case 'psthbinsize'
            psthbinsize = varargin{option+1};
        case 'minthresh'                          % ripple SD minimum threshold
            minthresh = varargin{option+1};
        case 'maxvelocity' %max velocity at which to exclude ripples
            maxvelocity = varargin{option+1};
        case 'pseudo'
            pseudodelay = varargin{option+1};     % impose offset (in sec) for error trials, to mimic delay in error. e.g. Egypt had 500 ms delay in reward
        case 'input'
            inputtimes = varargin{option+1};
        case 'output'
            outputtimes = varargin{option+1};
        case 'doripples'
            doripples = varargin{option+1};
        otherwise
            error(['Option ''', varargin{option}, ''' not defined']);
    end
end

out.index = index;

% establish times vector (suited for histc, where center bin straddles time 0)
% FOR RASTER
times = (-window(1)-0.5*binsize):binsize:(window(2)+0.5*binsize); %for use with histc binedges
timescenters = times(1:(end-1))+binsize/2;

% FOR PSTH
times2 = (-window(1)-0.5*psthbinsize):psthbinsize:(window(2)+0.5*psthbinsize); %for use with histc binedges
times2centers = times2(1:(end-1))+psthbinsize/2;

out.times = timescenters;
out.psthtimes = times2centers;

pos_samp = 30.03; %sampling rate of position (frames per s)  %set this as fixed because even if it varies a little bit (sometimes 29.97); if it is very consistent, use: 1/(pos{index(1)}{index(2)}.data(2,1)-pos{index(1)}{index(2)}.data(1,1));
pos_win=round([window(1)*pos_samp window(2)*pos_samp]);
veltimes = -window(1):1/pos_samp:window(2);
out.veltimes = veltimes;

reward = rewardinfo{index(1)}{index(2)};
spiketimes = spikes{index(1)}{index(2)}{index(3)}{index(4)}.data(:,1);

% if exist(cellinfo{index(1)}{index(2)}{index(3)}{index(4)}.type) <<to use this, must add cellinfo as an function argument
% celltype = cellinfo{index(1)}{index(2)}{index(3)}{index(4)}.type;
% end
%cellarea = cellinfo{index(1)}{index(2)}{index(3)}{index(4)}.area;

% Exclude spike times that occur during CA1 ripples
if ~isempty(excludeperiods)
    for s = 1:length(spiketimes)
        for e = 1:size(excludeperiods,1)
            if excludeperiods(e,1)<spiketimes(s) && spiketimes(s)<excludeperiods(e,2)
                spiketimes(s)=NaN;
            end
        end
    end
end

out.spiketimes = spiketimes;
%collect error and reward triggers
% first gather unique reward triggers (sweep through all epochs)
wells = [0 1 2 3 4 5];  % can manually specify to insure all wells are included even if not visited, or set to empty to run the loop
%             for d=index(1) %1:20
%                 for e=1:30
%                     try
%                         wells = [wells ; rewardinfo{d}{e}(:,1)];
%                     catch
%                     end
%                 end
%             end
%      wells = unique(wells(~isnan(wells)));

% construct trigmatrix --  error (0) and reward (1) triggers and inbound (10)-outbound (11)
% also initialize output
trigmatrix = [];
for w=1:length(wells)
    trigmatrix = [trigmatrix ; wells(w) 0 10];  % well error inbound
    trigmatrix = [trigmatrix ; wells(w) 1 10];  % well reward inbound
    trigmatrix = [trigmatrix ; wells(w) 0 11];  % well error outbound
    trigmatrix = [trigmatrix ; wells(w) 1 11];  % well reward outbound
end
%Note that some trigmatrix rows (such as correct inbound at a non-home well) are impossible trial types - those cells in the
%trigger times array will always be empty.  They are only there for parallel structure.

% initialize output matrix
out.trigmatrix = trigmatrix;
out.trigmatrixdescript = 'well #  //  error (0) or reward (1) // inbound (10) or outbound (11)';
if exist('celltype','var')
    out.type = celltype;
end
out.spikeraster = {};
out.spikepsth = {};
out.ripples = {};
out.rewardinputtimes = {};
out.rewardoutputtimes = {};
out.triggertimes = {};
out.rewardoutcomes = {};
out.velocity = {};
for t=1:size(trigmatrix,1)
    out.rewardinputtimes{t} = [];
    out.rewardoutputtimes{t} = [];
    out.triggertimes{t} = [];
    out.spikeraster{t} = [];
    out.spikepsth{t} = [];
    out.ripples{t} = [];
    out.rewardoutcomes{t} = [];
    out.velocity{t} = [];
end

if doripples
    %%%% This is bascially a recapitulation of getripples
    % for this epoch, construct array where each element represents whether
    % ripple is present on cell's parent tetrode (1 ms timesteps).
    % 	r = ripples{index(1)}{index(2)}{index(3)};
    % --- to use ripplescons ---
    r = eventcons{index(1)}{index(2)}{1};  % 1 is for CA1
    timevec = r.timerange(1):0.001:r.timerange(end);
    nrip = zeros(size(timevec));
    % 	    % apply the minthresh threhsold
    rvalid = find(r.maxthresh > minthresh);
    % ---- for pulling ripples from individual CA1 tets (old) ----
    % [rippletimes, ripplestdout] = ms_getripples([index(1) index(2)], ripplescons, cellinfo,tetinfo,'tetfilter','(isequal($area, ''dCA1''))', 'minstd',minthresh,'minrip',2);
    % ----
    % 	    rtimes = [r.starttime(rvalid) r.endtime(rvalid)];
    % create another parallel vector with bordering times for zeros
    tmp_rippletimes = [r.starttime(rvalid) r.endtime(rvalid)];
    %find ripples where velocity>maxvelocity and eliminate those
    for ripno = 1:size(tmp_rippletimes,1)
        if pos{index(1)}{index(2)}.data(lookup(tmp_rippletimes(ripno,1),pos{index(1)}{index(2)}.data(:,1)),5)>maxvelocity
            tmp_rippletimes(ripno,:) = NaN;
        end
    end
    nan = find(isnan(tmp_rippletimes(:,1)));
    tmp_rippletimes(nan,:)=[];
    rippletimes = tmp_rippletimes;
    nrtimes = [(rippletimes(:,1) - 0.00001) (rippletimes(:,2) + 0.00001)];
    rtimes = reshape(rippletimes', length(rippletimes(:)), 1);
    rtimes(:,2) = 1;
    nrtimes = [r.timerange(1) ; reshape(nrtimes', ...
        length(nrtimes(:)), 1) ; r.timerange(2)];
    nrtimes(:,2) = 0;
    % create a new list with all of the times in it
    tlist = sortrows([rtimes ; nrtimes]);
    % use interp to create a set of ones and zeros for each time
    % and add to nrip to get a cumulative count of the number of
    % ripples per timestep
    nrip = interp1(tlist(:,1), tlist(:,2), timevec, 'nearest');
    ripstamps = nrip .* timevec;
    %     [rippletimes, ripplestdout] = getripples([index(1) index(2)], ripples, cellinfo,'cellfilter','(isequal($area, ''CA1'') && ($numspikes > 20))', 'minstd',minthresh,'minrip',2);
else
    ripstamps = [];
end

% collect trigger times
% iterate over each reward trigger, iterating by trigger type
for t=1:size(trigmatrix,1)
    
    % if pseudo specified, then check if error trial
    %   - if so, institute the known offset (e.g. 0.5 sec)
    if trigmatrix(t,2)==0
        offset = pseudodelay;
    elseif trigmatrix(t,2)==1
        offset = 0;
    end
    
    row = rowfind(trigmatrix(t,:),reward(:,[1 3 4]));
    while row ~= 0
        % transcribe reward time in .rewardtimes
        out.rewardinputtimes{t} = [out.rewardinputtimes{t} ; reward(row,5)];
        out.rewardoutputtimes{t} = [out.rewardoutputtimes{t} ; (reward(row,2)+offset*10000)];
        % transcribe reward outcome in .rewardoutcomes
        out.rewardoutcomes{t} = [out.rewardoutcomes{t} ; reward(row,3)];
        % strike off trigger
        reward(row,[1 3 4])=[NaN NaN NaN];
        % look for another
        row = rowfind(trigmatrix(t,:),reward(:,[1 3 4]));
    end
    
end

triggers=[];

%determine whether trigger is surrounded by enough time <maxvelocity,
%up until end of window - if not, exclude it

for t=1:size(trigmatrix,1)
    vel=[];
    veltimes = [];
    % Define triggering events as reward input or output times
    if inputtimes==1
        triggers = out.rewardinputtimes{t}/10000;
        reward_col = 5;
    elseif outputtimes==1
        triggers = out.rewardoutputtimes{t}/10000;
        reward_col = 2;
    end
    
    if ~isempty(triggers)  %remove triggers too close to the beginning or end
        while (triggers(1)-pos{index(1)}{index(2)}.data(1,1))<window(1)
            triggers(1) = [];
        end
        while ~isempty(triggers) && (pos{index(1)}{index(2)}.data(end,1)-triggers(end)) < window(2)
            triggers(end) = [];
        end
    end
    
    if ~isempty(triggers)
        for ii=1:length(triggers)
            currtime = triggers(ii);
            % collect velocity corresponding to the window
            velind = lookup(currtime,pos{index(1)}{index(2)}.data(:,1));
            vel_start=velind-pos_win(1);  % index of window start in pos vector
            vel_end=velind+pos_win(2);    % index of window end in pos vector
            if vel_start>=1 && vel_end<=length(pos{index(1)}{index(2)}.data(:,1))
                tmpvel =pos{index(1)}{index(2)}.data(vel_start:vel_end,5);  %9th column = smoothed vel, 5th column = unsmoothed vel
                tmptimes = (pos{index(1)}{index(2)}.data(vel_start:vel_end,1)-currtime); %Set new pos times to 0 at reward
            end
            %lookup the velocities 1.5 second after trigger through the end
            %(should be <4 if rat actually came to a stop)
            testperiod = tmpvel(lookup(1.5,tmptimes):end);
            
            %note that 30 position indices corresponds to 1 s - exclude
            %trials where >=0.5s of time in the window was spent at
            %velocities>maxvelocity
            rewardind = lookup(triggers(ii)*10000,reward(:,5));
            
            if sum(testperiod>4)>=15    %tmpvel(end)>(maxvelocity+1)
                triggers(ii) = NaN;
                continue
                %also only include trials where there is at least sum(window) sec
                %from the last trigger to this one
            elseif rewardind~=1 &&  (reward(rewardind,reward_col)-reward(rewardind-1,reward_col))/10000 < (window(1)+window(2))
                triggers(ii) = NaN;
                disp(['trial excluded for ' num2str(index)])
                continue
            else
                vel =[vel; tmpvel'];  %9th column = smoothed vel, 5th column = unsmoothed vel
            end
        end
        triggers(isnan(triggers))=[];
    end
    
    %transcribe trigger times after exclusions
    out.triggertimes{t} = [out.triggertimes{t} triggers];
    
    % Now collect spike times around remaining triggers
    if ~isempty(triggers)
        for ii=1:length(triggers)
            
            % if using outputtimes and pseudo specified, then check if error trial
            %   - if so, institute the known offset (e.g. 0.5 sec)
            if trigmatrix(t,2)==0 && outputtimes==1
                offset = pseudodelay;
            elseif trigmatrix(t,2)==1 && outputtimes==1
                offset = 0;
            end
            
            % collect spike and ripple times
            if inputtimes==1
                spikeraster = histc(spiketimes,triggers(ii) + times);
                spikepsth = histc(spiketimes,triggers(ii) + times2);
                riphist = histc(ripstamps,triggers(ii) + times);
            elseif outputtimes==1
                
                spikeraster = histc(spiketimes,triggers(ii) + offset + times);
                spikepsth = histc(spiketimes,triggers(ii) + offset + times2);
                riphist = histc(ripstamps,triggers(ii) + offset + times);
            end
            out.spikeraster{t} = [out.spikeraster{t} ; spikeraster(1:end-1)']; % the end of the hist is always a 0, so chop it off
            out.spikepsth{t} = [out.spikepsth{t} ; spikepsth(1:end-1)'];
            out.ripples{t} = [out.ripples{t} ; riphist(1:end-1)];
        end
    end
    out.velocity{t} = [out.velocity{t} ; vel];
    
end

end











