
%% To do
%make this compatibl with multi-antimal F output.. currently animals has to
%be a single animal
%% Dashboard
close all
runFilterFramework = 1;
; saveFilterOutput = runFilterFramework;
; loadFilterOutput = 0;
EpochMean = 0;
; resaveFilterOutput = 0;
plotfigs = 0;
; plotSpec_EpochMean = 0;
; plotSpec_allEpochs = 0;
; savefigs = 0;
; pausefigs = 1;
%% ---------------- plotting params --------------------------
colorSet = 'DR1';
win = [.5 .5]; %in seconds
indwin = win*1500;
%% ---------------- Data Filters --------------------------
filtfunction = 'riptrigspectrogram';
animals = {'JZ1'};
days = [1:6];
eventtype = 'rippleskons';
eventSourceArea = 'ca1';
% epochEnvironment = 'openfield'; %wtrack, wtrackrotated, openfield, sleep
epochType = {'sleep', 'run', 'run'};
epochEnvironment = {'sleep', 'wtrack', 'openfield'};% 'wtrack'; %wtrack, wtrackrotated, openfield, sleep
eventSourceArea = 'ca1'; %ca1, mec, por
ripAreas = {'ca1', 'mec', 'por'};
ntAreas = {'ca1', 'sub', 'mec', 'por', 'v2l', 'ref'};

consensus_numtets = 1;   % minimum # of tets for consensus event detection
minstdthresh = 3;        % STD. how big your ripples are
exclusion_dur = 0.5;  % seconds within which consecutive events are eliminated / ignored
minvelocity = 0;
maxvelocity = 4;
%% ---------------- Paths and Title strings ---------------------------------------------------
investInfo = animaldef(lower('Demetris'));
filtOutputDirectory = sprintf('%s%s/', investInfo{2}, filtfunction);
figdirectory = sprintf('%s%s/', investInfo{4}, filtfunction);
filenamesave = sprintf('%s%sSTD%d_%s_%s_%s_D%s', eventSourceArea, eventtype, minstdthresh, strjoin(epochEnvironment,'-'), cell2mat(animals), cell2mat(LFPtypes),strjoin(arrayfun(@(x) num2str(x),days,'un',0),'-'));
filename = sprintf('%s_%s.mat', filtfunction, filenamesave);
filename = sprintf('riptrigspec_%s%s_%s_%sTets_%s.mat', eventSourceArea, eventtype, epochEnvironment, eventSourceArea, cell2mat(animals));
filenameTitle = strrep(filename,'_', '\_');
%% ---------------- Run FIlter ---------------------------------------------------
if runFilterFramework == 1;
    eptypeEnv = [epochType; epochEnvironment];
    epochfilter = sprintf('((isequal($type, ''%s'')) && (isequal($environment, ''%s''))) || ',eptypeEnv{:});
    yuck = strfind(epochfilter, '||');
    epochfilter = epochfilter(1:yuck(end)-1);  %cut off the trailing '||'
    %     iterator = 'multitetrodeanal'; %multitetrodeanal
    epochfilter =    sprintf('(isequal($type, ''%s'')) && (isequal($environment, ''%s''))',epochType, epochEnvironment); %'isequal($type, ''run'') && (isequal($environment, ''MultipleW''))'; %%'(isequal($type, ''sleep''))'; %%%&& isequal($descript, ''post-allruns''))';%   %%% %'isequal($type, ''run'') && isequal($environment, ''WTrackA'') && (($exposure>=1) && ($exposure<=10))';  %
    iterator = 'epocheeganal';
    %     tetfilter = sprintf('(isequal($area,''%s''))', eventSourceArea);
    tetfilter = sprintf('(isequal($area,''%s'')) || ', ntAreas{:});
    yuck = strfind(tetfilter, '||');
    tetfilter = tetfilter(1:yuck(end)-1); %cut off the trailing '||'
    timefilter{1} = {'getconstimes', '($cons == 1)',[eventSourceArea,eventtype],1,'consensus_numtets',consensus_numtets,...
        'minstdthresh',minstdthresh,'exclusion_dur',exclusion_dur,'minvelocity',minvelocity,'maxvelocity',maxvelocity};
    %----------F = createfilter('animal', animals, 'days', dayfilter,'epochs', epochfilter, 'excludetime', timefilter, 'eegtetrodes',tetfilter,'iterator', iterator);--------
    F = createfilter('animal', animals, 'days', days,'epochs', epochfilter, 'excludetime', timefilter, 'eegtetrodes',tetfilter,'iterator', iterator);
    %----------f = setfilteriterator(f, funcname, loadvariabl   es, options)--------
    F = setfilterfunction(F, ['dfa_' filtfunction], {'eeg', [eventSourceArea,eventtype]},'eventtype',eventtype);
    tic
    F = runfilter(F);
    F.filterTimer = toc;
    F.worldDateTime = clock;
    F.dataFilter = struct('animal', animals, 'days', days,'epochs', epochfilter, 'excludetime', timefilter, 'eegtetrodes',tetfilter,'iterator', iterator);
end
%% ---------------- Save Filter Output ---------------------------------------------------
if saveFilterOutput == 1;
    if ~isdir(filtOutputDirectory);
        mkdir(filtOutputDirectory);
    end
    save(sprintf('%s/%s',filtOutputDirectory, filename), 'F','-v7.3');
    disp(sprintf('filteroutput saved to %s/%s',filtOutputDirectory, filename))
end
%% ---------------- Load Filter Output ---------------------------------------------------
if loadFilterOutput == 1;
    load(sprintf('%s/%s',filtOutputDirectory, filename))
    disp(sprintf('filteroutput loaded: %s/%s',filtOutputDirectory, filename))
end
%% ---------------- Collect across tetrodes ---------------------------------------------------
if EpochMean == 1;
    %collect all indices into dayeptet
    indices = [];
    for i=1:length(F.output{1})  % iterate over epochs
        indices = [indices ; F.output{1}(i).index];
    end
    dayeptet = unique(indices(:,:),'rows');
    tets = unique(indices(:,3));
    
    %Collect data on individual tetrodes
    F.tetout = struct;
    for t = 1:length(tets)
        inds = dayeptet(dayeptet(:,3)==tets(t),:);
        F.tetout(t).indices = inds;
        F.tetout(t).t = F.output{1}(1).t;
        F.tetout(t).f = F.output{1}(1).f;
        %initialize
        F.tetout(t).S=[];
        F.tetout(t).numtrigs=[];
        for c = 1:length(F.output{1})
            if rowfind(F.output{1}(c).index,inds)~=0
                F.tetout(t).S = cat(3,F.tetout(t).S,F.output{1}(c).S);
                F.tetout(t).numtrigs = [F.tetout(t).numtrigs; numel(F.output{1}(1).triggers)];
            end
        end
    end
    %compute the average per tetrode
    F.meanspecs = cell(1,length(tets));
    for t=1:length(tets)
        F.meanspecs{t} = nanmean(F.tetout(t).S,3);
    end
end

% add tags from tetinfo struct as a field in the F.tetout
%load the tetinfostruct
animalinfo = animaldef(animals{1});
load([animalinfo{2},animalinfo{1} 'tetinfo']);
supareainfo = cellfetch(tetinfo, 'suparea');
areainfo = cellfetch(tetinfo, 'area');
subareainfo = cellfetch(tetinfo, 'subarea');

%scavange the tags from each tetrode
for ci = 1:length(F.tetout)
    tmpInds = F.tetout(ci).indices(1,:); %just grab the first Day Ep Tet ind... each tetrode should only be in one area for the duration of the exp
    
    targetInd = find(ismember(supareainfo.index, tmpInds, 'rows'));
    F.tetout(ci).suptag = supareainfo.values{targetInd};
    
    targetInd = find(ismember(areainfo.index, tmpInds, 'rows'));
    F.tetout(ci).areatag = areainfo.values{targetInd};
    
    targetInd = find(ismember(subareainfo.index, tmpInds, 'rows'));
    F.tetout(ci).subtag = subareainfo.values{targetInd};
end

%% ---------------- reSave Filter Output ---------------------------------------------------
if resaveFilterOutput == 1;
    if ~isdir(sprintf('%s/filter_output/riptrigspectrogram/', filtOutputDirectory));
        mkdir(sprintf('%s/filter_output/riptrigspectrogram/', filtOutputDirectory));
    end
    %save the entire workspace for filter provenance
    save(sprintf('%s/filter_output/riptrigspectrogram/riptrigspec_%s%s_%s_%sTets_%s.mat',filtOutputDirectory, eventSourceArea, eventtype, epochEnvironment, eventSourceArea, cell2mat(animals)), 'F','-v7.3');
end
%% plot
% if plotfigs

%% ---------------- Plot Averaged Across Epochs ---------------------------------------------------
if plotSpec_EpochMean == 1
    %     if ~exist('meanspecs');
    %         meanspecs = cell(1,length(tets));
    %         for t=1:length(tets)
    %             meanspecs{t} = mean(F.tetout(t).S,3);
    %         end
    %     end
    %the total height and length of the figure dictates the sizing of the subplots
    figLeft = 1;
    figBottom = 1;
    figHeight = 400;
    figlength = 1800;
    sfSpacing = .05;
    
    sfStartLeft = sfSpacing;
    sfEndRight = 1 - sfSpacing;
    sfAllLength = sfEndRight - sfStartLeft;
    
    sfAllBottom = .2;
    sfAllTop = .8 - sfAllBottom;
    
    figHandle = figure('Position', [figLeft, figBottom, figlength, figHeight+figBottom]);
    sfLength = sfAllLength/length(F.tetout); %normalize
    sfLeftBottom= linspace(sfStartLeft, sfEndRight-sfLength, length(F.tetout));
    %     subRightTop =linspace(subLength, 1, length(tets));
    %     subLength = floor(figlength/length(tets));
    %     subLeftBottom= floor(linspace(figLeft, figlength-subLength, length(tets)));
    %     subRightTop = floor(linspace(subLength, figlength, length(tets)));
    filenameTitle = sprintf('riptrigspec\\_%s%s\\_%s\\_%sTets\\_%s', eventSourceArea, eventtype, epochEnvironment, eventSourceArea, cell2mat(animals));
    %if we are dealing with cortical tetrodes, sort by layer and color code
    %the subfigure titles by layer
    if F.tetout(1).suptag == 'ctx';
        titletags = [];
        for t = 1:length(F.tetout)
            titletags = [titletags; F.tetout(t).subtag];
        end
        [titleColortags, sortInds]  = sort(titletags);
    else
        sortInds = 1:length(F.tetout);
        titletags = ones(1,length(F.tetout));
    end
    mycolors = lines;
    for t=1:length(F.tetout)
        %         subplot(1, length(tets),t, 'Position', [1, 1, 10, 10])
        %         subplot(1, length(tets),t,'Position', [subLeftBottom(t)+1, figBottom+1, subRightTop(t), figHeight])
        positionVector = [sfLeftBottom(t), sfAllBottom, sfLength, sfAllTop];
        subplot('Position', positionVector)
        %         subplot('positionVector', [subLeftBottom(t), 0, subRightTop(t), 1])
        imagesc(F.tetout(sortInds(t)).t,F.tetout(sortInds(t)).f,F.meanspecs{sortInds(t)}',[-0.2,2.5])
        set(gca,'YDir','normal')
        %         plot(rand(3,4))
        %         axis xy
        colormap('jet')
        %             xlabel('seconds from rip start','FontSize',12,'FontWeight','bold','Color','k')
        %             title({[sprintf('%s d%de%d %sRip(%d) Triggered %s-Spectrogram',animal, day, epoch, rippleregion, currrip, LFPtype{LFP})];[sprintf('timestamp: %16.f', centerripStartTime)]; [sprintf('peak nstds: %d', ripout{srcRegionind}{day}{epoch}.maxthresh(currrip))]},'FontSize',12,'FontWeight','bold')
        if t == 1
            xlab = xlabel('Time since start of detected SWR');
            set(xlab, 'horizontalAlignment', 'left')
            ylabel('Frequency')
            %             title(sprintf('tet %d', tets(t)))
            title(sprintf('nT:%d sA:%d', F.tetout(sortInds(t)).indices(1,3), titletags(sortInds(t))), 'Color',mycolors(titletags(sortInds(t)),:))
        else
            set(gca, 'YTick', []);
            set(gca,'XTick', [])
            set(gca, 'XTickLabel', [])
            %             title(sprintf('%d', tets(t)))
            %             title(sprintf('tet:%d lr:%d', F.tetout(sortInds(t)).indices(1,3), F.tetout(sortInds(t)).subtag), 'Color',mycolors(F.tetout(sortInds(t)).subtag,:))
            title(sprintf('nT:%d sA:%d', F.tetout(sortInds(t)).indices(1,3), titletags(sortInds(t))), 'Color',mycolors(titletags(sortInds(t)),:))
        end
        if t == length(F.tetout)
            %             colorbar
            supertitle(filenameTitle)
        end
    end
end
%% ---------------- Plot By Epochs ---------------------------------------------------
if plotSpec_allEpochs == 1;
    for e=1:3:length(F.output{1})
        figure
        imagesc(F.output{1}(e).t,F.output{1}(e).f,mean(F.output{1}(e).S,3)',[-0.2,2.5])
        axis xy
        xlabel('Time since start of detected SWR')
        ylabel('Frequency')
        colormap('jet')
        colorbar
    end
end

% end