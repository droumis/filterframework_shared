function [wellsdio, rewardinfo] = createrewardinfo_lintrack(directoryname,prefix,days,epochs, varargin)
%% This function parses DIO into a rewardinfo struct for the linear track.  In rewardinfo{day}{epoch}:
%       Column 1 = well visited
%       Column 2 = output timestamp (reward delivery) in NSpike time units 
%       Column 3 = trial logic (1 = correct, 0 = incorrect)
%       Column 4 = trajectory (inbound = 10, outbound = 11)
%       Column 5 = input timestamp (nosepoke)
% based on sj_findwellsfromdio1 [wellsdio] =
% 
%example use:  [wellsdio, rewardinfo] = createrewardinfo_lintrack('/opt/data40/mari/Bro','bro',[1 3 4 5],[2 4 6]);
%dio bit numbers not necessary if specified per animal below, but can specify as varargin: 
% (...,'inputdios',[well0 well1],'outputdios',[well0 well1]);

%% From DIO, gets well start and end for all unique trajectories
format long
lowercasethree = '';
%set variable options
for option = 1:2:length(varargin)-1
    switch varargin{option}
        case 'lowercasethree'
            lowercasethree = varargin{option+1};
        case 'inputdios'
            inputdios = varargin{option+1};
        case 'outputdios'
            outputdios = varargin{option+1};
    end
end

for day=days,
    
    dsz = '';
    if (day < 10)
        dsz = '0';
    end
    
    DIOfile = sprintf('%s/%sDIO%02d.mat', directoryname, prefix, day);
    load(DIOfile);

    if exist('inputdios','var') && exist('outputdios','var')
        % if DIO input and output bits specified in varargin, use those
        usedios = inputdios;
        outdios = outputdios;
    else
    % IF you know order of DIOs - enter here. List as well 1, well 2
    
    %% EXAMPLE
    %usedios = [24,23];      % INPUTS: well 1, well 2
    %outdios=[8,7];            % OUTPUTS: well 1, well 2
    

    end
    
    for epoch=epochs,
        
        %Initialize
        well_start=[]; well_end=[]; well_startend=[]; wellseq_curr=[]; welltrigtime_curr=[]; tmpwells=[];
        
        DIO1=DIO{day}{epoch}{usedios(1)}; nones = length(DIO1.pulsetimes); % Well 0: Bits 24 in
        DIO2=DIO{day}{epoch}{usedios(2)}; nzeros = length(DIO2.pulsetimes); % Well 1: Bits 25 in
      
        % Output wells
        DIOout1=DIO{day}{epoch}{outdios(1)}; nones_out = length(DIOout1.pulsetimes); % Well 0: Bits 8 out
        DIOout2=DIO{day}{epoch}{outdios(2)}; nzeros_out = length(DIOout2.pulsetimes); % Well 1: Bits 9 out

        out_trigt = [DIOout1.pulsetimes(:,1);DIOout2.pulsetimes(:,1)];
        [out_trigtx, out_ord] = sort(out_trigt);
        out_wells = [1*ones(nones_out,1);0*ones(nzeros_out,1)];
        out_rewwells = out_wells(out_ord);
        
        wells = [1*ones(nones,1);0*ones(nzeros,1)];
        trigt = [DIO1.pulsetimes(:,1);DIO2.pulsetimes(:,1)];
        [tfx,ti] = sort(trigt);
        wfx=wells(ti);
        swit=find(diff(wfx)~=0)+1;
        well_start=[wfx(1);wfx(swit(1:end-1))];
        well_end=[wfx(swit)];
        well_startend=[well_start,well_end];
        % For saving
        wellsdio{day}{epoch}=well_startend;
        
        wellseq_curr=[wfx(1);wfx(swit)];
        welltrigtime_curr=[tfx(1);tfx(swit)];
        
        % IF First Well is 1, then first outbound was captured. Keep it aside temporarily
         % Remove first outbound if it exists. Put back later
         tmpwells = well_startend; flag_firstout=0;
         if wellseq_curr(1)==1, 
             tmpwells(1,:)=[]; 
             flag_firstout=1;
         end
%         if wellseq_curr(1)~=1,
%             tmpseq=[1;tmpseq];
%             first = [1, well_startend(1,1)];
%             tmpwells = [first; tmpwells];
%         end
            
        %% Find outbound (rat starts at well 8 (1) on linear track)
        outbound_stidx = find(tmpwells(:,1)==1);
        outbound_wellstend = tmpwells(outbound_stidx,:);
        outbound_logic = zeros(length(outbound_wellstend),1);
        corridx = find( (tmpwells(outbound_stidx,2)~=tmpwells(outbound_stidx,1))...
            & (tmpwells(outbound_stidx,2)~=1) );                                    %Finds trips where end well does not equal start well
        correct_out = outbound_stidx(corridx);
        wrong_out = setdiff(outbound_stidx,correct_out);
        %% Find inbound (rat starts at well 9 (0) on linear track)
        inbound_stidx = find(tmpwells(:,1)~=1);
        inbound_wellstend = tmpwells(inbound_stidx,:);
        % Find out which inbound is correct, and return the logic
        inbound_logic = zeros(length(inbound_wellstend),1);
        corridx = find( (tmpwells(inbound_stidx,2)~=tmpwells(inbound_stidx,1))...
            & (tmpwells(inbound_stidx,2)==1) );
        correct_in = inbound_stidx(corridx);
        wrong_in = setdiff(inbound_stidx,correct_in);
        %inbound_logic(corr) = 1;
        
        % IMPORTANT - Put first Outbound Back In If It Was Removed - Have to Push all indexes By 1
        %if wellseq_curr(1)~=1,
        if wellseq_curr(1)==1,
            first = well_startend(1,:);
            tmpwells = [first; tmpwells];           
            % Push out and in idxs by 1
            correct_out = correct_out+1; correct_out=[1;correct_out];
            correct_in = correct_in+1;
            wrong_out = wrong_out+1;
            wrong_in = wrong_in+1;
            outbound_stidx = outbound_stidx+1; outbound_stidx=[1;outbound_stidx]; 
            inbound_stidx = inbound_stidx+1;
        end
                
        all_correct = sort([correct_out; correct_in]);
        all_wrong = sort([wrong_out; wrong_in]);
        allwell_curr = sort([all_correct;all_wrong]);
        rewtime_curr = welltrigtime_curr(all_correct+1);
        norewtime_curr = welltrigtime_curr(all_wrong+1);
        rewardedwell_seq = well_startend(all_correct,2); unrewardedwell_seq = well_startend(all_wrong,2);
        wellseq_reconstructed = wellseq_curr(allwell_curr+1);
        
        % Now prepare outputs
        wellseq_out =  wellseq_curr(2:end); 
        welltrigtime_out = welltrigtime_curr(2:end);
        logic = zeros(size(wellseq_out)); logic(all_correct)=1;
        
        trajseq_out = 10*ones(size(wellseq_out)); %10=inbound
        trajseq_out(outbound_stidx) = 11;
        
        
        % Compare rewarded time with those from Output DIOS
        output_idxs = lookup(rewtime_curr, out_trigtx); x=out_trigtx(output_idxs);
        difft = (rewtime_curr - out_trigtx(output_idxs))./10; % in msec; /10000 for sec
        replace = find(abs(difft)<2100); % Difference is usually on the orders of tens of ms
        rewtime_curr(replace) = x(replace);
        % Replace in the long sequence
        a=find(logic==1);
        a_replace=a(replace);
        welltrigtime_out_updated = welltrigtime_out;
        welltrigtime_out(a_replace) = x(replace);
        
        
        rewardinfo_curr = [wellseq_out, welltrigtime_out, logic, trajseq_out, welltrigtime_out_updated]; 
        rewardinfo{day}{epoch} = rewardinfo_curr;
        
    end
      % Columns of rewardinfo:
        % 1 = end well of trajectory (column 2 of wellsdio)
        % 2 = output time; reward delivery time on correct trials, input time on incorrect trials
        %       NOTE: for trials with (e.g.) 500 ms delay between input and output, make sure abs(difft)>500, to 
        %       include output times that are within something greater than 500 ms from the input time
        % 3 = logic; 1 is correct, 0 is incorrect
        % 4 = trial type; 10 is inbound, 11 is outbound
        % 5 = input (nosepoke) time on all trials
    
    
    
    
    eval([lowercasethree,'wellsdio = wellsdio;']); eval([lowercasethree,'rewardinfo = rewardinfo;']);
    if (directoryname(end) ~= '/')
        animdirect = [directoryname '/'];
    end
    eval(['save ',animdirect,prefix,'rewardinfo',dsz,num2str(day),' ',lowercasethree,'rewardinfo']);
    
end


