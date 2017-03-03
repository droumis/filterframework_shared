function [wellsdio, rewardinfo] = createrewardinfo_multW_trackbacks(directoryname,prefix,days,epochs,rewarddelay,varargin)
%% This function parses DIO into a rewardinfo struct for the multiple W track.  In rewardinfo{day}{epoch}:
%       Column 1 = well visited
%       Column 2 = output timestamp (reward delivery) in NSpike time units 
%       Column 3 = trial logic (1 = correct, 0 = incorrect)
%       Column 4 = trajectory type (inbound = 10, outbound = 11)
%       Column 5 = input timestamp (nosepoke)
% Written by Mari Sosa, inspired by sj_findwellsfromdio1_wtrack written by Shantanu Jadhav.

%example uses:  
    % To get rewardinfo for the rewarded sequences: [wellsdio, rewardinfo] = createrewardinfo_multW('/opt/data40/mari/Eli','eli',18,[3 5 7],2,'sequencetype',['S2';'S1';'S2']);
    % To get rewardinfo for the opposite/unrewarded sequences: [wellsdio, rewardinfo] = createrewardinfo_multW('/opt/data40/mari/Eli','eli',18,[3 5 7],2,'sequencetype',['S1';'S2';'S1'],'oppositesequence',1);
%dio bit numbers not necessary if specified per animal below, but can specify as varargin: 
% (...,'inputdios',[0 1 2 3 4 5],'outputdios',[left center right]); (left, center, right of reality, not pos reconstruct, e.g. 2 3 4 -- can only specify as varargin if doing 1 sequence at a time)

% rewarddelay = delay in seconds from nosepoke to reward. 
    % e.g. 1, if 1 second for all wells
	% e.g. [1,2,3] if wells vary and delays are either 1, 2, or 3 seconds

%% From DIO, gets well start and end for all unique trajectories
format long
%set variable options
for option = 1:2:length(varargin)-1
    switch varargin{option}
        case 'inputdios'
            inputdios = varargin{option+1};
        case 'outputdios'
            outputdios = varargin{option+1};
        case 'sequencetype'
            sequencetype = varargin{option+1};
        case 'oppositesequence'
            oppositesequence = varargin{option+1};
    end
end

if ~exist('oppositesequence','var')
    oppositesequence = 0; %default if unspecified
end

for day=days,
    
    dsz = '';
    if (day < 10)
        dsz = '0';
    end
    
    DIOfile = sprintf('%s/%sDIO%02d.mat', directoryname, prefix, day);
    posfile = sprintf('%s/%spos%02d.mat', directoryname, prefix, day);
    load(DIOfile);
    load(posfile);
    
        
    if ~exist('sequencetype','var')
        sequencetype = ['S1'; 'S1'; 'S1']; %default if unspecified, applies to S1 acquisition days
    end
    if exist('inputdios','var')
        % if DIO input bits specified in varargin, use those
        indios = inputdios;   
    else
        %If you know order of INPUT DIOs, enter here: list DIO cells in order 0, 1, 2, 3, 4, 5 (where 2 is center on S2, 3 is center on S1)

        % EXAMPLE
        switch prefix
            case 'fab'
                % FABIO
                indios = [32, 31, 30, 29, 28, 27];
                animflag = 4;
        end
    end
    
    % Define outputdios by sequence type
    for runnum = 1:length(epochs)
        epoch = epochs(runnum);
        if strcmp(sequencetype(runnum,:),'S1')
            L = 2; C = 3; R = 4;
            if exist('outputdios','var')
                outdios = outputdios;
            else % EXAMPLE
                if animflag == 4
                outdios = [8,7,6];       % for S1, Fab       % OUTPUTS: list DIO cells in order left, center, right (of reality, not pos reconstruct) (e.g. 2 3 4) 
                end
            end
        elseif strcmp(sequencetype(runnum,:),'S2')
             L = 1; C = 2; R = 3;
             if exist('outputdios','var')
                 outdios = outputdios;
             else % EXAMPLE
                 if animflag == 4
                     outdios = [5,8,7];       % for S2, Fab
                 end
             end


     % Output wells
        %First check if any output DIOs are empty
        filled_outdios = outdios;
        for i = 1:length(outdios)
            if isempty(DIO{day}{epoch}{outdios(i)})  % for empty DIO cells
                errormessage = sprintf('DIO %d is empty',outdios(i));
                disp(errormessage)
                % Exclude the empty DIO from the DIO input list
                filled_outdios(i) = NaN;
            end
        end
        % Delete NaNs
        filled_outdios(isnan(filled_outdios))=[];
        
        % Collect output trigger times and well IDs for non-empty output DIOs
        wellID = zeros(1,length(filled_outdios));
        out_trigtimes = [];
        out_wells = [];
        for j = 1:length(filled_outdios)
%             DIOin{j} = DIO{day}{epoch}{filled_indios(j)};
            if filled_outdios(j) == outdios(1)
            wellID(j) = L;
            elseif filled_outdios(j) == outdios(2)
                wellID(j) = C;
            elseif filled_outdios(j) == outdios(3)
                wellID(j) = R;
            end
            out_trigtimes = [out_trigtimes; DIO{day}{epoch}{filled_outdios(j)}.pulsetimes(:,1)]; 
        out_wells = [out_wells; wellID(j)*ones(size(DIO{day}{epoch}{filled_outdios(j)}.pulsetimes,1),1)];
        end 
        
        out_all = [out_trigtimes out_wells]; 
        % Sort well visits chronologically
        out_all = sortrows(out_all,1);
        
     % Input wells
        % First check if any are empty - sj_diodayprocess returns the cell
        % as empty if there are no timestamps for that well for the entire
        % day.  If the empty cell is at the end of the DIO array, it will
        % not exist.
        filled_indios = indios;
        for i = 1:length(indios)
            if indios(i)>length(DIO{day}{epoch}) % for nonexistent DIO cells
                errormessage = sprintf('DIO %d does not exist',indios(i));
                disp(errormessage)
                % Exclude the nonexistent DIO from the DIO input list
                filled_indios(i) = NaN;
                continue
            elseif isempty(DIO{day}{epoch}{indios(i)})  % for empty DIO cells
                errormessage = sprintf('DIO %d is empty',indios(i));
                disp(errormessage)
                % Exclude the empty DIO from the DIO input list
                filled_indios(i) = NaN;
            end
        end
        % Delete NaNs
        filled_indios(isnan(filled_indios))=[];
        
        % Collect input trigger times and well IDs for non-empty DIOs
        DIOin = cell(1,length(filled_indios));
        wellID = zeros(1,length(filled_indios));
        in_trigtimes = [];
        in_wells = [];
        for j = 1:length(filled_indios)
            DIOin{j} = DIO{day}{epoch}{filled_indios(j)};
            if filled_indios(j) == indios(1)
            wellID(j) = 0;
            elseif filled_indios(j) == indios(2)
                wellID(j) = 1;
            elseif filled_indios(j) == indios(3)
                wellID(j) = 2;
            elseif filled_indios(j) == indios(4)
                wellID(j) = 3;
            elseif filled_indios(j) == indios(5)
                wellID(j) = 4;
            else
                wellID(j) = 5;
            end
            in_trigtimes = [in_trigtimes; DIOin{j}.pulsetimes(:,1)]; 
        in_wells = [in_wells; wellID(j)*ones(size(DIOin{j}.pulsetimes,1),1)];
        end
        
        in_all = [in_trigtimes in_wells]; 
        % Sort well visits chronologically
        in_all = sortrows(in_all,1);
        
        %% Remove erroneous input triggers, BUT KEEP TRACKBACK ERRORS 
        % 1. To use as a check later, find whether a trigger is the first
        % one at that well: 0 indicates that it is NOT the first input
        first_input = [1; diff(in_all(:,2))]; % pad with a 1 at the beginning, since 1st recorded input of epoch has to be the first input at that well 
        % 2. Look up all input trigger times in position data 
        posind=lookup(in_all(:,1)/10000,pos{day}{epoch}.data(:,1));
        % 3. Find the max Y distance traveled between the each trigger
        % and its previous trigger - assumes position is oriented with
        % track arms along the Y axis
        maxposchange = NaN; %pad with NaN for the first recorded input of the epoch
            for i = 2:length(posind)
                trajectory = pos{day}{epoch}.data(posind(i-1):posind(i),3); % all Y position points between triggers
                currYpos = pos{day}{epoch}.data(posind(i),3); % Y pos of the current trigger
                poschange = currYpos-trajectory; % subtract to get the distance traveled at each point
                maxposchange = [maxposchange; max(abs(poschange))]; % find the maximum distance traveled on the trajectory between triggers
            end
        %sanity checks of maxposchange against well IDs and first_input
        check1 = [maxposchange in_all(:,2)];
        check2 = [maxposchange first_input];
        % 4. To include trackback errors: keep input triggers where rat traveled >= 48 in the Y direction
        % between triggers. Cutoff = 48 = halfway down the arm, since that's when rat's tail clears
        %the well and there is no more possibility of an erroneous trigger 
        keep_inputs = [in_all(1,:); in_all(maxposchange>=48,:)];
        keep_maxposchange = [NaN; maxposchange(maxposchange>=48)];
        % Create a matrix of starts and ends of each trajectory, in which 3rd column has input timestamps of trajectory end well  
        well_startend=[keep_inputs(1:end-1,2) keep_inputs(2:end,2) keep_inputs(2:end,1)]; 
        % Re-Add first captured "inbound" trial (center well visit that starts the session)
        well_startend = [NaN keep_inputs(1,2) keep_inputs(1,1);well_startend];
        % For saving
        wellsdio{day}{epoch}=well_startend;
        
        %% Mark well visits as correct or incorrect according to the rewarded sequence
        
        %Find outbound trials (should always start with the first trial captured)
        centervisits = find(well_startend(:,1)==C); %finds trials where start well is the center well
        first_outbound = well_startend(centervisits(1),1:2);
        all_outbound = well_startend(centervisits,:);
        % Set up logic vector for outbound trials
        outbound_logic = zeros(size(all_outbound,1),1);
        % Assign 1st outbound as correct if either L or R arm is visited
        if first_outbound(:,2)==L || first_outbound(:,2)==R
            outbound_logic(1)=1;
        end
        % Assign remaining outbounds as correct if they alternate L and R
        correct_outbound = find(all_outbound(2:end,2)~=all_outbound(1:end-1,2)); % first find alternations (current end well is different than previous end well)
        correct_outbound = correct_outbound+1; % Push indices by 1 to align to all_outbound, since we exclude the 1st trial
        outbound_logic(correct_outbound)=1;  % assign alternations as correct
        nonrewarded = wellID((wellID~=L & wellID~=R));  % now find alternations that end at wells that are not L and R
        % Revert logic to 0 for "correct" alternations to wells that are not L or R
        incorrect_index = find(ismember(all_outbound(correct_outbound,2),nonrewarded)); 
        revert = correct_outbound(incorrect_index);
        outbound_logic(revert)=0;
        %put it all together with input timestamps and outbound label "11"
        OUTS = [all_outbound outbound_logic 11*ones(length(outbound_logic),1)]; 
        
        % NOTE: the above section will include as correct any outbound
        % visits to L or R (home-adjacent) wells where the previous inbound trial was initiated at a nonrewarded well
        %i.e. if a non-sequence well is visited on an outbound trial, on the NEXT outbound trial either home-adjacent arm is correct
        
        %Find inbound trials 
        notcenter = find(well_startend(:,1)~=C); %finds trials where start well is anything but center well
        all_inbound = well_startend(notcenter,:);
        % Set up logic vector for inbound trials
        inbound_logic = zeros(size(all_inbound,1),1);
        % Assign inbounds as correct if they end in the center well 
        correct_inbound = find(all_inbound(:,2)==C);
        inbound_logic(correct_inbound)=1;
        %put it all together with input timestamps and inbound label "10"
        INS = [all_inbound inbound_logic 10*ones(length(inbound_logic),1)]; 
        
        % Combine and sort inbounds and outbounds
        all_trials = [INS; OUTS];
        all_trials = sortrows(all_trials,3);
        
        %% Special case updates to logic
        % (A) Assigns as correct any outbound home-adjacent visits where
        % prev prev was anything but the current well (accounts for prev
        % prev inbound errors to nonrewarded wells), unless this outbound
        % follows a trackback error.
        
        % (B) Also assigns as INCORRECT any outbound home-adjacent visits where prev prev well was the same as the current well, even if this was on an incorrect inbound.
                
        %(C)  After the above are completed, find outbound trials that follow
        % trackback errors, and scan back in time to assign logic based on
        % most recent non-trackback. (trackbacks on inbound trials should
        % already be assigned as incorrect, since they don't end in the
        % home arm).
       
        find_L = find(all_trials(:,2)==L);
        find_R = find(all_trials(:,2)==R);
        for p = 1:length(find_L)
            prev_ind = find_L(p)-1;
            prevprev_ind = find_L(p)-2;
            if prev_ind ~=0 && prevprev_ind~=0
                prev_well = all_trials(prev_ind,2);
                prevprev_well = all_trials(prevprev_ind,2);
                % (A) update corrects that are not following a trackback
                if prev_well==C && prevprev_well~=L && prevprev_well~=C
                    if all_trials(find_L(p),4) == 0
                    disp(['update index' ' ' num2str(find_L(p))]); %displays indices of all_trials that actually updated
                    end
                    all_trials(find_L(p),4)=1;
                % (B) update incorrects
                elseif prev_well==C && prevprev_well==L
                    if all_trials(find_L(p),4) == 1
                    disp(['update index' ' ' num2str(find_L(p))]); %displays indices of all_trials that actually updated
                    end
                    all_trials(find_L(p),4)=0;
                % (C) update incorrect trials after trackbacks 
                % - correct trials should already be labeled as correct
                elseif prev_well==C && prevprev_well==C %this may not work if there is a legitimate trackback within the first few trials
                    scanwell = all_trials(prev_ind,2); %start with prev_well, C
                    while scanwell==C && prev_ind~=0
                        prev_ind = prev_ind-1;
                        if prev_ind==0
                            break
                        end
                        scanwell = all_trials(prev_ind,2);
                    end
                    if scanwell==L
                        if all_trials(find_L(p),4) == 1
                            disp(['update index' ' ' num2str(find_L(p))]); %displays indices of all_trials that actually updated
                        end
                    all_trials(find_L(p),4)=0;
                    end
                end
            end
        end
        for p = 1:length(find_R)
            prev_ind = find_R(p)-1;
            prevprev_ind = find_R(p)-2;
            if prev_ind ~=0 && prevprev_ind~=0
                prev_well = all_trials(prev_ind,2);
                prevprev_well = all_trials(prevprev_ind,2);
                % (A) update corrects that are not following a trackback
                if prev_well==C && prevprev_well~=R && prevprev_well~=C
                    if all_trials(find_R(p),4) == 0
                    disp(['update index' ' ' num2str(find_R(p))]); %displays indices of all_trials that actually updated
                    end
                    all_trials(find_R(p),4)=1;
                % (B) update incorrects
                elseif prev_well==C && prevprev_well==R
                    if all_trials(find_R(p),4) == 1
                    disp(['update index' ' ' num2str(find_R(p))]); %displays indices of all_trials that actually updated
                    end
                    all_trials(find_R(p),4)=0;
                % (C) update incorrect trials after trackbacks 
                % - correct trials should already be labeled as correct
                elseif prev_well==C && prevprev_well==C %this may not work if there is a legitimate trackback within the first few trials
                    scanwell = all_trials(prev_ind,2); %start with prev_well, C
                    while scanwell==C && prev_ind~=1
                        prev_ind = prev_ind-1;
                        scanwell = all_trials(prev_ind,2);
                    end
                    if scanwell==R
                        if all_trials(find_R(p),4) == 1
                            disp(['update index' ' ' num2str(find_R(p))]); %displays indices of all_trials that actually updated
                        end
                    all_trials(find_R(p),4)=0;
                    end
                end
            end
        end
        
        
        
        
        %% Add output timestamps
        
        % First add column for potential outputs
        input_times = all_trials(:,3);
        all_trials = [all_trials input_times];
        % Compare input timestamps with output DIO timestamps
        
        output_inds = lookup(input_times, out_all(:,1)); 
        possible_outputs = out_all(output_inds,1);
        % Find difference of inputs and possible outputs in seconds
        difft = (possible_outputs-input_times)./10000; 
        % Find output times that fall within the max reward delay,accounting for additional ~20 ms
        replace = find(abs(difft)<(max(rewarddelay)+0.1)); 
        % Replace input timestamps with output timestamps on correct trials in all_trials 
            % (don't need to specify correct trials - should already
            % align by looking up output times, which have to be correct)
        all_trials(replace,6) = possible_outputs(replace,1);
        % Sanity check - can display if desired
        checkreplace = [all_trials(:,6)-all_trials(:,3) all_trials(:,4)];
        checkall = [keep_maxposchange all_trials];
        
        %% Restructure all_trials into existing rewardinfo structure from Shantanu
        
        % Columns:
        % 1 = end well of trajectory
        % 2 = output (reward delivery) time on correct trials
        % 3 = logic; 1 is correct, 0 is incorrect
        % 4 = trial type; 10 is inbound, 11 is outbound
        % 5 = input (nosepoke) time on all trials
        rewardinfo_curr = [all_trials(:,2) all_trials(:,6) all_trials(:,4) all_trials(:,5) all_trials(:,3)];
        rewardinfo{day}{epoch} = rewardinfo_curr;
        
      end  
        % Save
        if (directoryname(end) ~= '/')
        animdirect = [directoryname '/'];
        end
        
        if oppositesequence %rename saved variable as oppseqrewardinfo, and filename as oppseqrewardinfo
            oppseqrewardinfo = rewardinfo;
            filename = sprintf('%s/%soppseqrewardinfo%d%d.mat',directoryname,prefix,str2num(dsz),day);
             save(filename,'oppseqrewardinfo');
        else %save rewardinfo for the rewarded sequence
    filename = sprintf('%s/%srewardinfo%d%d.mat',directoryname,prefix,str2num(dsz),day);
filename2 = sprintf('%s/%swellsdio%d%d.mat',directoryname,prefix,str2num(dsz),day);
    save(filename,'rewardinfo');
    save(filename2,'wellsdio');
        end
    
   
end