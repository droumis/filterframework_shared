%% Plot event (such as ripple) -triggered spike rasters and psths, and calculate significance of ripple modulation for each cell
% MS updated Feb 2017; originally adapted from Kenny's dfskk_eventtriggeredspiking

%%%%% GENERAL PARAMETERS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
eventconsname = 'vca1ripplescons';
anim = 'Ger';               % capitalized prefix
animals = {'Geronimo'};     % animal(s)
epochtype = 'run';

runscript = 1; 
usedayfilter = 1;
dayfilter = 5;              % enter days to filter here
postprocessing = 1;             
calcsig = 1;                % calculate significance of each cell's SWR-modulation

plot_singlerasters = 1;     % plot individual cell rasters and psths
savefigs = 0;               % save each figure
savedatastruct_tofile = 0;  % save post-processed data to file

%%% SPARSIFY %%%
% Option to randomly delete a certain proportion of spikes for each cell
% to match a control firing rate, during significance calculation
sparsify = 0;

%%%%% RUNSCRIPT processing parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if runscript
    animals_torun = animals; 
    disp(eventconsname)
    window = [0.5 0.5];         % size of psth window (in sec)
    binsize = .001;             % size of bins (in sec)
    frbinsize = 0.02;
    time = (-window(1)-0.5*binsize):binsize:(window(2)+0.5*binsize);
end
%%%%% POST-processing parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if postprocessing
    animals_toprocess =  animals;  %can change if only want to process 1 animal at a time
end 

%%%%% consensus timefilter parameters %%%%%
switch eventconsname
    case 'dca1ripplescons'
        datadir = sprintf('/opt/data40/mari/%s/',anim);     % the datadir to store processed files
        TF = 1; %'(isequal($validripple, 1))';              % tetfilter used in analysis function to transcribe postmatrix -- ripple detection has already occurred in ripplescons
        consensus_numtets = 3;                              % minimum # of tets for consensus event detection
        exclude_ripples = 0;                                % obviously 0; might change if event = wave gamma
        welldist = [];
        minthresh = 3;                                      % how big your ripples/gammaf are
        exclusion_dur = 0.5;                                % seconds within which consecutive events are eliminated / ignored
        maxvelocity = 4;                                    % max head velocity to accept events
        minvelocity = 0;
        timefilter = { {'getconstimes', '($cons == 1)', 'dca1ripplescons',1,...
            'consensus_numtets',consensus_numtets,'minthresh',minthresh,'exclusion_dur',exclusion_dur,'maxvelocity',maxvelocity} };
    case 'vsubripplescons' 
        datadir = sprintf('/opt/data40/mari/%s/',anim);                 
        TF = 1; %'(isequal($validripple, 1))';              
        consensus_numtets = 1;                              
        exclude_ripples = 0;  
        welldist = [];
        minthresh = 3;                                      
        exclusion_dur = 0.5;                                
        maxvelocity = 4;                                    
        minvelocity = 0;
        timefilter = { {'getconstimes', '($cons == 1)', 'vsubripplescons',1,...
            'consensus_numtets',consensus_numtets,'minthresh',minthresh,'exclusion_dur',exclusion_dur,'maxvelocity',maxvelocity} };
        case 'vca3ripplescons' 
        datadir = sprintf('/opt/data40/mari/%s/',anim);        
        TF = 1; %'(isequal($validripple, 1))';       
        consensus_numtets = 1;   
        exclude_ripples = 0;  
        welldist = [];
        minthresh = 3;          
        exclusion_dur = 0.5;             
        maxvelocity = 4;
        minvelocity = 0;
        timefilter = { {'getconstimes', '($cons == 1)', 'vca3ripplescons',1,...
            'consensus_numtets',consensus_numtets,'minthresh',minthresh,'exclusion_dur',exclusion_dur,'maxvelocity',maxvelocity} };
case 'vca1ripplescons' 
        datadir = sprintf('/opt/data40/mari/%s/',anim);        
        TF = 1; %'(isequal($validripple, 1))';       
        consensus_numtets = 1;   
        exclude_ripples = 0;  
        welldist = [];
        minthresh = 3;         
        exclusion_dur = 0.5;             
        maxvelocity = 4;
        minvelocity = 0;
        timefilter = { {'getconstimes', '($cons == 1)', 'vca1ripplescons',1,...
            'consensus_numtets',consensus_numtets,'minthresh',minthresh,'exclusion_dur',exclusion_dur,'maxvelocity',maxvelocity} };
end
      

 %% Run DF Script    

 if runscript
     
     % Cell filter
     cellfilter = '(isequal($area, ''NAc'') && ($numspikes > 100))';  
     
     % Iterator
     iterator = 'singlecellanal';
     
     for animal = animals_torun
         
         animalinfo = animaldef(animal{1});
         animaldir = animalinfo{2};
         animalprefix = animalinfo{3};
         if strcmp(epochtype,'sleep')
             epochfilter =  '(isequal($type, ''sleep''))'; 
         elseif strcmp(epochtype,'run')
             epochfilter =  '(isequal($type,''run'') && isequal($environment,''MultipleW''))';
         end
         
         % Filter Creation
         if usedayfilter
             f = createfilter('animal', animal,'days',dayfilter,'epochs', epochfilter, 'cells', cellfilter, 'excludetime', timefilter, 'iterator', iterator);
         else
             % ca1f = createfilter('animal', animal, 'epochs', epochfilter, 'cells', ca1cellfilter, 'excludetime', timefilter, 'iterator', iterator);
             f = createfilter('animal', animal,'epochs', epochfilter, 'cells', cellfilter, 'excludetime', timefilter, 'iterator', iterator);
         end
         
         % Set Analysis Function
         % ca1f = setfilterfunction(ca1f, 'dfakk_geteventtriggeredspiking', {'spikes',eventconsname,'pos','task','sleep','cellinfo'},animaldir,animalprefix,'TF',TF,'window',window,'binsize',binsize,'minthresh',minthresh,'maxvelocity',maxvelocity,'minvelocity', minvelocity,'consensus_numtets',consensus_numtets,'welldist',welldist);
         f = setfilterfunction(f, 'dfams_geteventtrigspiking', {'spikes',eventconsname,'pos','task','sleep','cellinfo'},animaldir,animalprefix,'TF',TF,'window',window,'binsize',binsize,'frbinsize',frbinsize,'minthresh',minthresh,'maxvelocity',maxvelocity,'minvelocity', minvelocity,'consensus_numtets',consensus_numtets,'welldist',welldist);
         
         % Run Analysis
         % ca1f = runfilter(ca1f);
         f = runfilter(f);
         
         % append runscript parameters to the raw output structs
         runscript_params = paramsstruct(eventconsname,timefilter,window,binsize,time,consensus_numtets,minthresh,exclusion_dur);
         %     ca1f.runscript_params = runscript_params;
         f.runscript_params = runscript_params;
         
         % if running for multiple regions and you want to save the output for future use
         % cd(datadir)
         % save(sprintf('Super%s_raw_%s_%s',eventconsname,animal{1}(1:3),date),'ca1f','ca2f','ca3f','dgf','nacf','animal','-v7.3');
         
     end
     
 end

 % f = nacf;

%% POSTPROCESSING

if postprocessing
   
    for animal = animals_toprocess
         
        % load cellinfo  (used below to gather number of clustered epochs)
        animalinfo = animaldef(animal{1});
        cellinfo = loaddatastruct(animalinfo{2}, animalinfo{3}, 'cellinfo');
        
        % iterate over regions
        for reg = 1 %1:6
            
            % load data if saved previously
%             regionscript  % to obtain rawvar
            cd(datadir)
%             filename = dir(sprintf('Super%s_raw_%s*',eventconsname,animal{1}(1:3)));  % if previously saved and want to load...
%             load([datadir filename.name],rawvar)
%             regionscript  % a script of "if" statements with color assignments, etc for various brain regions; see Kenny's
            
            % find unique cells (day-tet-cell)
            detc = [];
            if ~isempty(f.output)
                for ii=1:length(f.output{1})  % iterate over epochs
                    detc = [detc ;  f.output{1}(ii).index];
                end
                dtc = unique(detc(:,[1 3 4]),'rows');
            else
                continue
            end
            
            % identify the total # of epochs for each day
                % (used to initialize matrices properly below)
            day_numepochs = [];
            for d = unique(dtc(:,1))'
                day_numepochs = [day_numepochs ; d  length(cellinfo{d})];   %  [ day  <total # of epochs>]
            end
            
            % initialize output
            A = struct;
            
            % iterate through the animal's cells
            for ind=1:size(dtc,1)
                disp([animal{1}(1:3) ' ' num2str(dtc(ind,:))])
                %identify total number of epochs for this day
                numepochs = day_numepochs(day_numepochs(:,1)==dtc(ind,1),2);
                % initialize output .fields
                A(ind).dtc = dtc(ind,:);           
%                 A(ind).neighbor_principals = [];    % detc of principal units clustered at the same tetrode
                A(ind).epoch_types = [];            % indices here correspond to epoch #s
                A(ind).epoch_envs = [];             % indices here correspond to epoch #s
                % all detected events, regardless of epoch or state
                A(ind).psth = [];
                A(ind).frhist = [];
                A(ind).instantFR = [];
                A(ind).psthsum = [];
                A(ind).instantFRmean = [];
                A(ind).posteventmatrix = [];
                A(ind).eventduration = [];
                A(ind).eventtags = [];
                A(ind).eventtags_descript = '[ <epoch #>  <time of event>  < sleepc-occurring or not>]';
                A(ind).nospikes = 0;
                A(ind).noevents = 0;                  % number of events reported for the epochs in which the unit was clustered ("events experienced")
                A(ind).epochs = [];   
                % run epochs
                 A(ind).run_epochs = [];
                 A(ind).run_nospikes = nan(1,numepochs);
                 A(ind).run_noevents = nan(1,numepochs);
                % sleep epochs 
                 A(ind).sleep_epochs = [];
                 A(ind).sleep_nospikes = nan(1,numepochs);
                 A(ind).sleep_noevents = nan(1,numepochs); 
                % time vector for psth (bin centers)
                A(ind).time = [];
                 
                  % Consolidate data over each epoch
                for c=1:length(f.output{1})
                    % if day tet cell matches
                    if rowfind(dtc(ind,:),f.output{1}(c).index([1 3 4]))
                        % some older animals have ghost (not .type labelled) epochs -- ignore these
                        if ~strcmp(f.output{1}(c).epoch_type,'run') && ~strcmp(f.output{1}(c).epoch_type,'sleep')
                            disp('the .type field for an epoch in the task struct must be either "run" or "sleep" -- ignoring epoch')
                            continue
                        end
                        epochnum = f.output{1}(c).index(2);                        
                        % both-epoch outputs
                        A(ind).time = f.output{1}(c).time;   % (all time vectors are the same..)
                        A(ind).frtime = f.output{1}(c).frtime;
                        A(ind).psth = [A(ind).psth ; f.output{1}(c).psth];
                        A(ind).frhist = [A(ind).frhist; f.output{1}(c).frhist];
                        A(ind).instantFR = [A(ind).instantFR ; f.output{1}(c).instantFR];
                        A(ind).posteventmatrix = [A(ind).posteventmatrix ; f.output{1}(c).posteventmatrix];
                        A(ind).eventduration = [A(ind).eventduration ; f.output{1}(c).eventduration];
                        A(ind).eventtags = [A(ind).eventtags ; f.output{1}(c).eventtags];
                        A(ind).nospikes = A(ind).nospikes + f.output{1}(c).nospikes;
                        A(ind).noevents = A(ind).noevents + f.output{1}(c).noevents;
                        A(ind).epochs = [A(ind).epochs epochnum];
                        % epoch-by-epoch outputs
                        A(ind).epoch_types{epochnum} = f.output{1}(c).epoch_type;
                        A(ind).epoch_envs{epochnum} = f.output{1}(c).epoch_environment;                           
                        A(ind).epoch_nospikes(epochnum) = sum(f.output{1}(c).nospikes);
                        A(ind).epoch_noevents(epochnum) = sum(f.output{1}(c).noevents);
                        A(ind).epoch_noeventspikes(epochnum) = sum(sum(f.output{1}(c).psth));
                        % run vs. sleep outputs
                        if strcmp(f.output{1}(c).epoch_type,'run')
                            A(ind).run_epochs = [A(ind).run_epochs  f.output{1}(c).eventtags(1,1)];
                            A(ind).run_nospikes(epochnum) = f.output{1}(c).nospikes;
                            A(ind).run_noevents(epochnum)  = f.output{1}(c).noevents;
                        elseif strcmp(f.output{1}(c).epoch_type,'sleep')
                            A(ind).sleep_epochs = [A(ind).sleep_epochs  f.output{1}(c).eventtags(1,1)];
                            A(ind).sleep_nospikes(epochnum) = sum(f.output{1}(c).nospikes);
                            A(ind).sleep_noevents(epochnum) = sum(f.output{1}(c).noevents);
                        
                        end
                    end
                end
                
                    % if a no data for this cell, ignore
                if isempty(A(ind).psth)
                    continue
                end
                
                % Convert .psth field to sparse matrix for smaller saved variable
                A(ind).psth = sparse(A(ind).psth);   
                %Sum psth
                A(ind).psthsum = sum(A(ind).psth,1);
                A(ind).instantFRmean = mean(A(ind).instantFR,1);
            end
        end
    end
end

                

%% RASTER PLOT %%%%%%%%%%%%%%%%%%%%%%%%%%

% Plotting options
mark_confirmed_sleep_ripples = 0;  % only use if "sleep" data structure was passed into the DFA
timecourse_colormap = 0;
plot_variance = 0;

region = 'NAc';

for c = 1:length(A)
    data = A(c);
    
    % Create smoothed PSTH of mean firing rate
    smoothing_length = 10;   % std of gaussian (in ms) used to smooth rip psth
    smoothing_width = round(smoothing_length*.001/f.runscript_params.binsize);   % smoothing width in number of bins
    kernel = gaussian(smoothing_width,smoothing_width*8);
    smoothedpsth = smoothvect(sum(full(data.psth),1)./(binsize*data.noevents),kernel);
    A(c).smoothedpsth = smoothedpsth;
    
    if plot_singlerasters
        
        h=figure;
        subplot(8,2,[1 3 5 7])
        hold on
        
        % before plotting the spike rasters themselves, overlay colored rectangles that indicate epochs
        alleps = unique(data.eventtags(:,1));
        alleps(isnan(alleps)) = [];
        
        for epno = alleps'
            firsttrial = rowfind(epno,data.eventtags(:,1));
            lasttrial = length(data.eventtags(:,1)) - rowfind(epno,flipud(data.eventtags(:,1))) + 1;
            
            % draw rectangle indicating RUN epoch
            if ismember(epno,data.run_epochs)
                patch([          data.time(1)*1000                  data.time(end)*1000          data.time(end)*1000     data.time(1)*1000],...
                    [           firsttrial - 0.5                    firsttrial - 0.5             lasttrial + 0.5        lasttrial + 0.5],...
                    [1  1  1]) %,'facealpha',.95)
            else
                % draw rectangle indicating REST epoch
                patch([          data.time(1)*1000                  data.time(end)*1000          data.time(end)*1000     data.time(1)*1000],...
                    [           firsttrial - 0.5                    firsttrial - 0.5             lasttrial + 0.5        lasttrial + 0.5],...
                    [.75 .75 .75]) %,'facealpha',.95)
            end
        end
        
        % (Optional) mark ripples that occur during periods of confirmed sleep
        if mark_confirmed_sleep_ripples
            % draw 0 line for sleepc ripples
            plot(zeros(size(data.psth,1),1),1:size(data.psth,1),'Color',[1 .4 .4],'linewidth',4)
            % draw patches corresponding to sleepc periods
            
            sleepc_eventvec = data.eventtags(:,3);
            list = vec2list(sleepc_eventvec,1:length(sleepc_eventvec));
            
            for pp = 1:size(list,1)
                % draw rectangle
                patch([          data.time(1)*1000                  data.time(end)*1000          data.time(end)*1000     data.time(1)*1000],...
                    [           list(pp,1) - 0.5                    list(pp,1) - 0.5             list(pp,2) + 0.5        list(pp,2) + 0.5],...
                    [.75 .75 .75],'edgecolor','none')  % ,'facealpha',.3
            end
        end
        
        % PLOT SPIKE RASTERS %%%%%%%%%%%%
        axisfontsize = 16;
        manual_clr = [0 0 0]; % set marker color if not specified by regionscript
        subplot(8,2,[1 3 5 7])
        if isempty(manual_clr)
            plotraster5(data.time,full(data.psth),85,clr,'burstisi',6,'postripplematrix',data.posteventmatrix)
        else
            plotraster5(data.time,full(data.psth),85,manual_clr,'burstisi',6,'postripplematrix',data.posteventmatrix)
        end
        
        % Plot epoch start and end lines (a bit redundant with rectangles,
        % but could change colors if desired)
        for epno = alleps'
            firsttrial = rowfind(epno,data.eventtags(:,1));
            lasttrial = length(data.eventtags(:,1)) - rowfind(epno,flipud(data.eventtags(:,1))) + 1;
            plot([data.time(1) data.time(end)],[firsttrial firsttrial],'-k','LineWidth',10)
            plot(data.time,lasttrial,'-k','LineWidth',10)
        end
        
        set(gca,'fontsize',axisfontsize,'fontweight','normal')
        
        % (Optional) overlay colormap of ripples over course of day
        if timecourse_colormap
            xposition = window(2)/2;
            ypositions = 1:size(data.psth,1);
            % create appropriate length colormap
            colormap_scaled = [];
            colormapvals = colormap;
            for chan = 1:3
                colormap_scaled(:,chan) = interp(colormapvals(:,chan),ceil(length(ypositions)/size(colormapvals,1)));
                colormap_scaled((colormap_scaled(:,chan) < 0),chan) = 0;
                colormap_scaled((colormap_scaled(:,chan) > 1),chan) = 1;
            end
            for q = 1:length(ypositions)
                plot(xposition,ypositions(q),'.','color',colormap_scaled(q,:));
            end
        end
        
        
        % SMOOTHED HISTOGRAM PLOT %%%%%%%%%%%%%%%%%%%%%%%%%
        subplot(8,2,[9 11])
        h = bar(data.time,smoothedpsth,'facecolor',[0 0 0]);
        
        % Set a background color if desired
        % bgcolor = [];
        % set(gca,'Color',bgcolor)
        
        axis tight
        hold on
        
        % x axis
        set(gca,'fontsize',axisfontsize,'fontweight','normal')
        xlabel('Time (sec)','fontsize',16,'fontweight','normal')
        ylabel('Firing Rate (Hz)','fontsize',16,'fontweight','normal')
        clear ylim;
        
        % limit
        ymaximum = max(sum(smoothedpsth,1));
        
        % plot 0 line
        if ymaximum < 10
            plot([0 0],[0 ymaximum+0.5],'linewidth',2,'Color',[1 .4 .4])
            ylim([0 ymaximum+0.2])
        elseif ymaximum > 10
            ymaximum2 = ceil(ymaximum)+2;
            plot([0 0],[0 ymaximum2],'linewidth',2,'Color',[1 .4 .4])
            ylim([0 ymaximum2])
        end
        
        % (Optional) plot variance - should look about the same as psth outline
        %smoothing_width = round(smoothing_length*.001/binsize);   % smoothing width in number of bins
        %kernel = gaussian(smoothing_width,smoothing_width*8);
        if plot_variance
            subplot(8,2,[14 15])
            psth_variance_smoothed = smoothvect(var(full(data.psth),1),kernel);
            var_mean = mean(psth_variance_smoothed);
            y_shift = ymaximum/2;
            y_scale = 0.5 * ymaximum/range(psth_variance_smoothed);
            psth_variance_toplot = (psth_variance_smoothed - var_mean) * y_scale + y_shift;
            h = plot(data.time,psth_variance_smoothed,'r','linewidth',3);
            %         set(gca,'Color',bgcolor)
            %         axis tight
        end
        
        titlestring=sprintf('%s %s %d %d %d (%s, %d)',...
            animal{1}(1:3),region,A(c).dtc,'riptrig-spikes', sum(full(data.psthsum)));
        
        [~,title_handle] = suplabel([titlestring],'t');
        set(title_handle,'Fontsize',18,'FontWeight','normal')
        set(gcf, 'renderer', 'zbuffer')
        if savefigs
            figfilename = sprintf('%sripplespikerasters/%s/%s_%d-%d-%d_%s.fig',datadir,epochtype,anim,A(c).dtc,eventconsname);
            saveas(h,figfilename)
        end
                close(gcf)
        
    end 
end

if savedatastruct_tofile
    cd(sprintf('%s/ripplemoddata',datadir))
    save(sprintf('%srippletrigspiking_%d_%s_100spikes.mat',animal{1}(1:3),day,eventconsname),'A','-v7.3');
end

%%  Calculate significance of modulation

close all

if calcsig
    
    % option to save ripplemod struct, below
    save_ripplemod = 0;
    daylist = unique(dtc(:,1))';
    tic
    nshuffles = 1000;
    % set up ripplemod struct for storing p values - currently saves all cells and days in one struct
    ripplemod = struct;
    ripplemod.days = daylist;
    ripplemod.event = eventconsname;
    ripplemod.celltype = cellfilter;
    ripplemod.epochtype = epochtype;
    
    % previously had an option to save this for each day, a la "eliripplemod04"
    % ^^ do we want this??
    %       for day = daylist
    %       tmpdtc =  dtc(dtc(:,1)==day,:);
    
    % Modvalues = [day, tet, cell, direction(pos=+1, neg=-1), pvalue, depth-of-modulation, mean-FR-in-window]
    ripplemod.modvalues = [dtc repmat(zeros(size(dtc,1),1),1,4)];
    
    for cellind = 1:size(dtc,1)      % (!!) if doing this for cherry-picked cells, enter cellind here
        %             if dtc(cellind,1) == day
        analyzing_cell = dtc(cellind,:);
        disp(['analyzing cell: ' num2str(analyzing_cell)])
        
        spikedata = full(A(cellind).psth);
        numspikes = sum(full(A(cellind).psthsum),2);
        
        % mean FR in the plotted window
        windowFR =  numspikes/((length(A(cellind).psthsum)*binsize)*size(A(cellind).psth,1));
        
        % (Optional) Randomly delete a certain proportion of spikes for each cell
        % to match a control firing rate, specified above.  "sparsify" must equal 1
        control = 0.5; %control FR to match
        if windowFR > control && sparsify
            spikeinds = find(spikedata);
            flip = rand(size(spikeinds));
            deletespikes = flip<=0.5;
            spikedata(spikeinds(deletespikes))=0;
            windowFR = sum(sum(spikedata,1),2)/((length(A(cellind).psthsum)*binsize)*size(A(cellind).psth,1))
            sparsified = 1;
            disp('sparsified')
        end
        
        % time range to analyze significance... could also use the average event duration, on either side of 0
        range = [0 0.2]; %mean(A(cellind).eventduration);
        % varRange = index for range(1) to range(2) in time series
        varRange = [lookup(range(1),A(cellind).time):lookup(range(2),A(cellind).time)];
        
        % set windowisrange to 1 to only shuffle spikes within the analysis range, not the whole -0.5 to 0.5
        windowisrange = 0;
        if windowisrange
            tmpspikedata=spikedata(:,varRange); %for shuffling only within the analysis range
            spikedata = tmpspikedata;
        end
        
        binwindow = size(spikedata,2);
        halfbinwindow = size(spikedata,2)/2;
        shufStep = [-halfbinwindow halfbinwindow]; % maximum time to shuffle spikes
        
        allShufPsth = zeros(nshuffles,binwindow);
        
        % smoothing should be the same as the plotted real data
        smoothing_length = 10;   % std of gaussian (in ms) used to smooth rip psth
        smoothing_width = round(smoothing_length*.001/f.runscript_params.binsize);   % smoothing width in number of bins
        kernel = gaussian(smoothing_width,smoothing_width*8);
        
        % circularly shuffle individual ripple-triggered spiketrains by r # of bins,
        for n = 1:nshuffles
            currshuffle = zeros(size(spikedata));
            for rip = 1:size(spikedata,1) % number of rips depends on number of clustered epochs for the cell
                %create an empty "spiketrain"
                dummyspikes = zeros(1,size(spikedata,2));
                %generate random value from the uniform distribution -0.5 to .5 (peri-ripple window in bins)
                r = round(shufStep(1) + (shufStep(2)-shufStep(1)).*rand(1,1));
                %find indices of spiketimes in current train
                spikeinds = find(spikedata(rip,:));
                shiftedspikeinds = spikeinds+r;
                % if times fall off the beginning or end, bring them around circularly
                circindend = find(shiftedspikeinds>size(spikedata,2));
                shiftedspikeinds(circindend)=shiftedspikeinds(circindend)-size(spikedata,2);
                circindstart = find(shiftedspikeinds<1);
                shiftedspikeinds(circindstart)=size(spikedata,2)-(1-shiftedspikeinds(circindstart));
                %fill in dummy spike train
                dummyspikes(shiftedspikeinds) = 1;
                currshuffle(rip,:) = dummyspikes;
            end
            % smoothed the shuffled psth
            smoothedcurrpsth = smoothvect(sum(currshuffle,1)./(binsize*size(spikedata,1)),kernel);
            %      figure
            %      h = bar(1:length(smoothedcurrpsth),smoothedcurrpsth,'facecolor',[0 0 0]);
            allShufPsth(n,:) = smoothedcurrpsth;
        end
        
        %take mean of shuffles
        meanShuf = mean(allShufPsth,1);
        
        %if shuffle was only done within the analysis range, use the whole thing
        if windowisrange
            meanShufRange = meanShuf;
        else
            % otherwise we analyze only the shuffled psth in "range"
            meanShufRange = meanShuf(varRange);
        end
        % figure
        % h = bar(1:length(meanShuf),meanShuf,'facecolor',[0 0 0]);
        
        %find summed squared distance of each shuffle
        meanvarShuf = zeros(nshuffles,1);
        for n=1:nshuffles
            if windowisrange
                allShufRange = allShufPsth(n,:);
            else
                allShufRange = allShufPsth(n,varRange);
            end
            meanvarShuf(n,:) = sum(((allShufRange-meanShufRange).^2),2);
        end
        
        %find summed squared distance of the real data psth from the mean of the shuffles, within the analysis range
        RealPsth = smoothvect(sum(spikedata,1)./(binsize*size(spikedata,1)),kernel);
        if windowisrange
            RealPsthRange = RealPsth;
        else
            RealPsthRange = RealPsth(varRange);
        end
        meanvarReal = sum((RealPsthRange-meanShufRange).^2);
        
        if sparsified
        figure
        h = bar(1:length(RealPsth),RealPsth,'facecolor',[0.5 0.5 0.5]);
        % figure
        % h = bar(varRange,RealPsthRange,'facecolor','b');
        % figure
        % h1 = bar(varRange,meanShufRange,'facecolor','k');
        end
        
        % p value is 1- (what fraction of shuffles have a variance that exceeds the real data)
        p = 1 - sum(meanvarShuf<meanvarReal)/nshuffles
        
        % direction of modulation
        if  sum(RealPsthRange-meanShufRange,2)<0
            direction = -1;
        else
            direction = 1;
        end
        
        ripplemod.modvalues(cellind,4) = direction;
        ripplemod.modvalues(cellind,5) = p;
        ripplemod.modvalues(cellind,6) = meanvarReal; %depth of modulation
        ripplemod.modvalues(cellind,7) = windowFR;
        
        %save ripplemod variable
        if save_ripplemod
            modname = sprintf('%sripplemoddata/%sripplemod_%s_%s_%s.mat',datadir,anim,epochtype,eventconsname,date);
            save(modname,'ripplemod','-mat')
        end
        
    end
    
    %         end
    ripplemod.modvalues
    %     end
    toc
end
