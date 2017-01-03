%%  Plots input- or output-DIO triggered spiking, where input is nosepoke and output is reward delivery. Either reward or error.
% inspired by dfskk_rewardtrigspiking
% relies on rewardinfo struct created by createrewardinfo_<typeoftrack>
% relies on "trigmatrix" output by dfams_getrewardtrigspiking

window = [4 4];
binsize = 0.001; %in sec, for raster
psthbinsize = 0.1; % in sec, for psth; larger bins for smoother psth

%Define trigger type
input=0;
output=1; % if 1, specify delay length as an argument to DFA to create "pseudo" output times for error trials

%Include ripples in raster?
dorip = 1;

%Plotting options
plotrewarderrorrasters = 1; %Plot spike rasters for reward vs. error trials
plotbywell = 0; %Plot rasters by well
plotpsth = 1; %Plot spike histograms
plotvelocity = 1; %Plot velocity with rasters
if plotvelocity == 1
    dovelocity = 1;
else
    dovelocity = 0;
end

datadir = '/opt/data40/mari/Eli/';

runscript = 1;
if runscript
    
    % Animal Selection
    animals = {'Eliot'};
    
    % Day Filter
    dayfilter = [19];
    
    % Epoch Filter
    %     epochfilter = 'isequal($type, ''run'') && (isequal($environment, ''WTrackA'') || isequal($environment, ''WTrackB''))';
    epochfilter = '(isequal($type, ''run'') && isequal($environment,''MultipleW''))';%(isequal($sequence, ''S1''))';
    
    % Time Filter
    timefilter{1} = {'kk_getconstimes','($cons==0)','dca1ripplescons',1,'consensus_numtets',2,'minthresh',3,'exclusion_dur',0,'maxvelocity',4,'minvelocity',0};
    timefilter{2} = {'kk_getconstimes','($cons==0)','vca3ripplescons',1,'consensus_numtets',1,'minthresh',3,'exclusion_dur',0,'maxvelocity',4,'minvelocity',0};
    % ^ exlcude times during consensus ripples
    
    %     ca1cellfilter = '(isequal($area, ''dCA1'') && ($numspikes > 20))'; %&& isequal($type, ''principal''))';    %% ($meanrate < 3))';
    %     ca2cellfilter = '(isequal($area, ''CA2'') && ($numspikes > 20))'; %&& isequal($type, ''principal''))';
    %     ca3cellfilter = '(isequal($area, ''vCA3'') && ($numspikes > 20))'; %&& isequal($type, ''principal''))';
    %     dgcellfilter = '(isequal($area, ''DG'') && ($numspikes > 20))'; %&& isequal($type, ''principal''))';
    naccellfilter = '(isequal($area, ''NAc'') && ($numspikes > 100))';
    
    % Iterator
    iterator = 'singlecellanal';
    
    %Number of regions to plot - set if iterating over multiple brain regions
    %     numregions = 3;
    
    % Filter Creation
    %     ca1f = createfilter('animal', animals, 'days', dayfilter,'epochs', epochfilter, 'cells', ca1cellfilter, 'excludetime', timefilter, 'iterator', iterator);
    %     ca2f = createfilter('animal', animals, 'days', dayfilter,'epochs', epochfilter, 'cells', ca2cellfilter, 'excludetime', timefilter, 'iterator', iterator);
    %     ca3f = createfilter('animal', animals, 'days', dayfilter,'epochs', epochfilter, 'cells', ca3cellfilter, 'excludetime', timefilter, 'iterator', iterator);
    %     dgf = createfilter('animal', animals, 'days', dayfilter,'epochs', epochfilter, 'cells', dgcellfilter, 'excludetime', timefilter, 'iterator', iterator);
    nacf = createfilter('animal', animals, 'days', dayfilter,'epochs', epochfilter, 'cells', naccellfilter, 'excludetime', timefilter, 'iterator', iterator);
    
    % Set Analysis Function
    %     ca1f = setfilterfunction(ca1f, 'dfams_getrewardtrigspiking', {'spikes','rewardinfo','ripplescons','pos','cellinfo','tetinfo'},'window',window,'binsize',0.001,'minthresh',2,'input',1,'prioritize_ripples',0,'dovelocity',1);
    %     ca2f = setfilterfunction(ca2f, 'dfams_getrewardtrigspiking', {'spikes','rewardinfo','ripples','pos','cellinfo','tetinfo'},'window',[5 5],'binsize',0.001,'minthresh',3,'input',1,'prioritize_ripples',0);
    %     ca3f = setfilterfunction(ca3f, 'dfams_getrewardtrigspiking', {'spikes','rewardinfo','ripplescons','pos','cellinfo','tetinfo'},'window',window,'binsize',binsize,'minthresh',3,'maxvelocity',4,'input',1,'doripples',0,'dovelocity',1);
    %     dgf = setfilterfunction(dgf, 'dfams_getrewardtrigspiking', {'spikes','rewardinfo','ripples','pos','cellinfo','tetinfo'},'window',[5 5],'binsize',0.001,'minthresh',3,'input',1,'prioritize_ripples',0);
    nacf = setfilterfunction(nacf, 'dfams_getrewardtrigspiking', {'spikes','rewardinfo','dca1ripplescons','pos'},'window',window,'binsize',binsize,'psthbinsize',psthbinsize,'minthresh',3,'maxvelocity',4,'input',input,'output',output,'doripples',dorip);
    
    % Run Analysis
    %     ca1f = runfilter(ca1f);
    %     ca2f = runfilter(ca2f);
    %     ca3f = runfilter(ca3f);
    %     dgf = runfilter(dgf);
    nacf = runfilter(nacf);
    
end





%% Collect across epochs

for reg=5 %1:numregions       % iterate over regions
    if reg==1
        region='dCA1';
        f=ca1f;
    elseif reg==2
        region='CA2';
        f=ca2f;
    elseif reg==3
        region='vCA3';
        f=ca3f;
    elseif reg==4
        region='DG';
        f=dgf;
    elseif reg==5
        region = 'NAc';
        f=nacf;
    end
    
    % Collapse cells across epochs within a day
    
    %collect all cell indices, ignoring specific epochs, into daytetcell
    indices = [];
    for i=1:length(f.output{1})  % iterate over epochs
        indices = [indices ; f.output{1}(i).index];
    end
    daytetcell = unique(indices(:,[1 3 4]),'rows');
    
    f.celloutput = struct;
    
    
    for ind=1:size(daytetcell,1)        % iterate over each unique cell
        
        % initialize entry
        f.celloutput(ind).daytetcell = daytetcell(ind,:);
        f.celloutput(ind).rewardinputtimes = cell(1,size(f.output{1}(1).trigmatrix,1));
        f.celloutput(ind).triggertimes = cell(1,size(f.output{1}(1).trigmatrix,1));
        f.celloutput(ind).rewardoutputtimes = cell(1,size(f.output{1}(1).trigmatrix,1));
        f.celloutput(ind).trigmatrix = f.output{1}(1).trigmatrix;
        f.celloutput(ind).trigmatrixdescript = f.output{1}(1).trigmatrixdescript;
        f.celloutput(ind).ripples = cell(1,size(f.output{1}(1).trigmatrix,1));
        f.celloutput(ind).spikeraster = cell(1,size(f.output{1}(1).trigmatrix,1));
        f.celloutput(ind).spikepsth = cell(1,size(f.output{1}(1).trigmatrix,1));
        f.celloutput(ind).rewardoutcomes = cell(1,size(f.output{1}(1).trigmatrix,1));
        f.celloutput(ind).velocity = cell(1,size(f.output{1}(1).trigmatrix,1));
        % initialize individual output matrices
        for m=1:length(f.celloutput(ind).spikeraster)
            f.celloutput(ind).spikeraster{m}=[];
            f.celloutput(ind).spikepsth{m}=[];
            f.celloutput(ind).ripples{m}=[];
            f.celloutput(ind).rewardinputtimes{m} = [];
            f.celloutput(ind).rewardoutputtimes{m} = [];
            f.celloutput(ind).triggertimes{m} = [];
            f.celloutput(ind).rewardoutcomes{m}=[];
            f.celloutput(ind).velocity{m}=[];
        end
        
        %collect data across epochs
        for c=1:size(indices,1)
            if daytetcell(ind,:)==indices(c,[1 3 4])
                f.celloutput(ind).times = f.output{1}(c).times;
                f.celloutput(ind).psthtimes = f.output{1}(c).psthtimes;
                f.celloutput(ind).veltimes = f.output{1}(c).veltimes;
                f.celloutput(ind).spiketimes = f.output{1}(c).spiketimes;
                %                 f.celloutput(ind).type = f.output{1}(c).type;
                for m=1:length(f.celloutput(ind).spikeraster)
                    % concatenate
                    f.celloutput(ind).spikeraster{m}=[f.celloutput(ind).spikeraster{m} ; f.output{1}(c).spikeraster{m}];
                    f.celloutput(ind).spikepsth{m}=[f.celloutput(ind).spikepsth{m} ; f.output{1}(c).spikepsth{m}];
                    f.celloutput(ind).rewardinputtimes{m} = [f.celloutput(ind).rewardinputtimes{m}; f.output{1}(c).rewardinputtimes{m}];
                    f.celloutput(ind).rewardoutputtimes{m} = [f.celloutput(ind).rewardoutputtimes{m}; f.output{1}(c).rewardoutputtimes{m}];
                    f.celloutput(ind).triggertimes{m} = [f.celloutput(ind).triggertimes{m}; f.output{1}(c).triggertimes{m}];
                    f.celloutput(ind).rewardoutcomes{m}=[f.celloutput(ind).rewardoutcomes{m} ; f.output{1}(c).rewardoutcomes{m}];
                    f.celloutput(ind).ripples{m}=[f.celloutput(ind).ripples{m} ; f.output{1}(c).ripples{m}];
                    f.celloutput(ind).velocity{m}=[f.celloutput(ind).velocity{m} ; f.output{1}(c).velocity{m}];
                end
            end
        end
    end
    
    super(reg) = f;
    
end

%% Save


% Filename = 'Egy_inputtrigspiking_allcells40.mat';
%
% save(Filename,'super', '-v7.3');
% disp('super saved!')

%% plot RASTERS: reward vs error trials

if plotrewarderrorrasters
    
    for reg=5
        
        if reg==1
            region='dCA1'; clr=[0 0 0];
        elseif reg==2
            region='CA2'; clr='k'; %[0 0.7 0];
        elseif reg==3
            region='vCA3'; clr = 'k'; %clr=[1 0 0];
        elseif reg==4
            region='DG';
            clr=[1 0 1];
        elseif reg==5
            region = 'NAc';
            clr =[0 0 0];
        end
        
        % plot cell-by-cell
        for c=1:length(super(reg).celloutput)
            
            data = super(reg).celloutput(c);
            daytetcell = data.daytetcell;
            A = figure;
            hold on
            
            for eventtype = [0 1]
                
                if eventtype == 1
                    flag = 'reward';
                    velplot = 3;
                    if plotvelocity
                        subplot(2,2,1)
                    else
                        subplot(1,2,1)
                    end
                    title({[region ' (' num2str(daytetcell) ')  ' eventtype ' rasters']; ...
                        flag},'fontsize',18,'fontweight','bold');
                elseif eventtype == 0
                    flag = 'error';
                    velplot = 4;
                    if plotvelocity
                        subplot(2,2,2) %subplot(1,2,2)
                    else
                        subplot(1,2,2)
                    end
                    title(flag,'fontsize',18,'fontweight','bold');
                end
                
                % collect all reward trial spike rasters (+ ripple 'rasters'), disregarding wells, in
                % chronological order
                taggedspikerasters = [];
                taggedriprasters = [];
                taggedvelocity = [];
                
                triggers=data.triggertimes;
                counter = 0;
                
                for trig = find(data.trigmatrix(:,2)==eventtype)'   % selected by reward or error
                    % rewardtime and raster
                    if ~isempty(triggers{trig})
                        taggedspikerasters = [ taggedspikerasters ; triggers{trig} data.spikeraster{trig} ] ;
                        if dorip
                            taggedriprasters = [ taggedriprasters ; triggers{trig} data.ripples{trig} ] ;
                        end
                        if plotvelocity
                            taggedvelocity = [ taggedvelocity ; triggers{trig} data.velocity{trig} ] ;
                        end
                        counter = counter + 1;
                    end
                end
                
                % exit if there were no plottable trials
                if counter == 0
                    continue
                end
                
                % sort into chronological order + untag
                taggedspikerasters = sortrows(taggedspikerasters,1);
                spikerasters = taggedspikerasters(:,2:end);
                if dorip
                    taggedriprasters = sortrows(taggedriprasters,1);
                    riprasters = taggedriprasters(:,2:end);
                end
                if plotvelocity
                    taggedvelocity = sortrows(taggedvelocity,1);
                    velocitytraces = taggedvelocity(:,2:end);
                    meanvelocity = mean(velocitytraces,1);
                    stdvelocity = std(velocitytraces, 1);
                end
                
                hold on
                % plot (no ripples)
                if dorip==0
                    msplotraster2(data.times,spikerasters,0,1,'Color',clr,'linewidth',2)
                    % plot (with ripples)
                elseif dorip==1
                    msplotraster3(data.times,spikerasters,riprasters,0,1,'Color',clr,'linewidth',2)
                end
                % plot well entry and reward lines
                if input==1
                    poke = 0;
                    rew = 2;
                elseif output==1
                    poke = -2;
                    rew = 0;
                end
                plot([poke poke],[0 size(spikerasters,1)],'Color',[0.5 0.5 0.5],'LineWidth',1)
                plot([rew rew],[0 size(spikerasters,1)],'Color',[0.9 0.7 0],'LineWidth',1)
                set(gca,'xtick',[-window(1):2:window(2)]);
                set(gca,'XLim',[-window(1) window(2)]);
                
                if plotvelocity
                    subplot(2,2,velplot)
                    hold on
                    plot(data.veltimes,meanvelocity,'k','LineWidth',2)
                    plot(data.veltimes,meanvelocity+stdvelocity,'Color',[0.7 0.7 0.7],'LineWidth',1)
                    plot(data.veltimes,meanvelocity-stdvelocity,'Color',[0.7 0.7 0.7],'LineWidth',1)
                    % plot well entry and reward lines
                    plot([poke poke],[0 50],'Color',[0.5 0.5 0.5],'LineWidth',1)
                    plot([rew rew],[0 50],'Color',[0.9 0.7 0],'LineWidth',1)
                    set(gca,'xtick',[-window(1):2:window(2)]);
                    set(gca,'XLim',[-window(1) window(2)]);
                    set(gca,'YLim',[0 50]);
                end
            end
            %             filename = sprintf('%srewardspikerasters/withvelocity/%s_NAc_%d-%d-%d_rewardspikeraster.fig',datadir,animals{1},daytetcell(1),daytetcell(2),daytetcell(3));
            %             saveas(A,filename)
        end
    end
end


%% Plot PSTHs: reward vs. error vs. all trials
if plotpsth
    % set up smoothing gaussian
    figsmoothing_length = 50;   % std of gaussian in ms
    smoothing_width = round(figsmoothing_length*.001/psthbinsize);   % smoothing width in number of bins
    figkernel = gaussian(smoothing_width,smoothing_width*8);
    
    allsmoothedmean = [];
    allsmoothedsem = [];
    % plot cell-by-cell
    for c= 1:length(f.celloutput)
        ymax = 0.5;
        data = f.celloutput(c);
        daytetcell = data.daytetcell;
        A = figure;
        hold on
        
        legdeets = [];
        for eventtype = [0 1 2] %0 error, 1 reward, 2 all
            
            if eventtype == 1
                flag = 'reward';
                lineclr = [0 0.7 0.4];
                if plotvelocity
                    subplot(2,1,1)
                else
                    subplot(1,1,1)
                end
                
            elseif eventtype == 0
                flag = 'error';
                lineclr = [0 0.3 0.6];
                if plotvelocity
                    subplot(2,1,1)
                else
                    subplot(1,1,1)
                end
                
            elseif eventtype ==2
                flag = 'all';
                lineclr = 'k';
                if plotvelocity
                    subplot(2,1,1)
                else
                    subplot(1,1,1)
                end
            end
            
            % collect all reward trial psth values, disregarding wells
            spikepsth = [];
            velpsth = [];
            counter = 0;
            
            if eventtype == 0 || eventtype == 1
                trigind = find(data.trigmatrix(:,2)==eventtype)';
            else
                trigind = 1:size(data.trigmatrix,1);
            end
            
            for trig = trigind   % selected by reward or error or all
                if ~isempty(data.triggertimes{trig})
                    spikepsth = [spikepsth; data.spikepsth{trig}];
                    if plotvelocity
                        velpsth = [velpsth; data.velocity{trig}];
                    end
                    counter = counter + 1;
                end
            end
            
            % exit if there were no plottable trials
            if counter == 0
                continue
            end
            
            % find mean and std of the velocity for current trial type
            meanvelpsth = mean(velpsth,1);
            stdvelpsth = std(velpsth, 1);
            
            %find mean and smooth for psth
            psthFR = spikepsth./psthbinsize;   % you need to do this trial by trial for the std dev
            meanpsth = mean(psthFR,1);
            sempsth = std(psthFR,0,1)/sqrt(size(psthFR,1)-1);
            
            %pad mean and sem so the ends don't "droop" once smoothed
            meanpad1 = repmat(meanpsth(1),1,10);
            meanpad2 = repmat(meanpsth(end),1,10);
            sempad1 = repmat(sempsth(1),1,10);
            sempad2 = repmat(sempsth(end),1,10);
            
            meanpsth = [meanpad1 meanpsth meanpad2];
            sempsth = [sempad1 sempsth sempad2];
            
            smoothedmean= smoothvect(meanpsth,figkernel);
            smoothedmean = smoothedmean(length(meanpad1)+1:(end-length(meanpad2)));
            smoothedsem = smoothvect(sempsth,figkernel);
            smoothedsem = smoothedsem(length(sempad1)+1:(end-length(sempad2)));
            
            allsmoothedmean = [allsmoothedmean; smoothedmean];
            allsmoothedsem = [allsmoothedsem; smoothedsem];
            
            % for printing figure legends
            legdeets{eventtype+1} = flag;
            legdeets{eventtype+4} = size(psthFR,1);
            
            % adjustable ymax for plotting
            if ~isempty(smoothedmean)
                newmax = ceil(max([smoothedmean+smoothedsem]));
                
            else
                newmax = 1;
            end
            if newmax>ymax
                ymax = newmax;
            end
            
            hold on
            % plot
            plot(data.psthtimes,smoothedmean,'Color',lineclr,'LineWidth',2)
            plot(data.psthtimes,smoothedmean+smoothedsem,'Color',lineclr,'LineWidth',1)
            plot(data.psthtimes,smoothedmean-smoothedsem,'Color',lineclr,'LineWidth',1)
            % plot well entry and reward lines
            if input==1
                poke = 0;
                rew = 2;
            elseif output==1
                poke = -2;
                rew = 0;
            end
            plot([poke poke],[0 max(smoothedmean+smoothedsem)+5],'Color',[0.5 0.5 0.5],'LineWidth',2)
            plot([rew rew],[0 max(smoothedmean+smoothedsem)+5],'--','Color',[0.9 0.7 0],'LineWidth',2)
            set(gca,'xtick',[-window(1):2:window(2)]);
            xlim([-window(1) window(2)]);
            ylim([0 ymax])
            ylabel('Firing rate (Hz)')
            
            if plotvelocity
                subplot(2,1,2)
                hold on
                plot(data.veltimes,meanvelpsth,'Color',lineclr,'LineWidth',2)
                plot(data.veltimes,meanvelpsth+stdvelpsth,'Color',lineclr,'LineWidth',1)
                plot(data.veltimes,meanvelpsth-stdvelpsth,'Color',lineclr,'LineWidth',1)
                set(gca,'xtick',[-window(1):2:window(2)]);
                set(gca,'XLim',[-window(1) window(2)]);
                set(gca,'YLim',[0 50]);
                plot([poke poke],[0 50],'Color',[0.5 0.5 0.5],'LineWidth',2)
                plot([rew rew],[0 50],'--','Color',[0.9 0.7 0],'LineWidth',2)
                ylabel('Velocity')
            end
        end
        legend(sprintf('%s/%s = %d/%d',legdeets{2},legdeets{1},legdeets{5},legdeets{4}),'Location','SouthOutside')
        H = suplabel(sprintf('%s %s [%s] Green=Reward, Blue=Error, Black=All; PSTH',animals{1},region,num2str(daytetcell)),'t');
        set(H,'FontSize',18,'FontWeight','bold')
        %             filename = sprintf('%srewarderror/%s_NAc_%d-%d-%d_rewarderrorpsth.fig',datadir,animals{1},daytetcell(1),daytetcell(2),daytetcell(3));
        %             saveas(A,filename)
    end
    
end

