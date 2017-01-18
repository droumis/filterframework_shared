
%---------------- ripple triggered spectrogram parameters--------------------------
% eventcons = 'ca1rippleskons';
eventtype = 'rippleskons';
eventarea = 'ca1';

consensus_numtets = 1;   % minimum # of tets for consensus event detection
minthresh = 2;        % STD. how big your ripples are
exclusion_dur = 0.5;  % seconds within which consecutive events are eliminated / ignored   
minvelocity = 0;
maxvelocity = 4; 

% ---------------- Data Filters ---------------------------------------------------
% Animal Selection
animals = {'JZ1'};
    
% Day Filter
dayfilter = 2;%:7; %6:16; %Fabio days 6-13 for S1 acquisition, 14-21 for switch
   
% Epoch Filter
epochfilter =    '(isequal($type, ''run'')) && (isequal($environment, ''wtrack''))'; %'isequal($type, ''run'') && (isequal($environment, ''MultipleW''))'; %%'(isequal($type, ''sleep''))'; %%%&& isequal($descript, ''post-allruns''))';%   %%% %'isequal($type, ''run'') && isequal($environment, ''WTrackA'') && (($exposure>=1) && ($exposure<=10))';  %

iterator = 'epocheeganal';

tetfilter = '(isequal($area,''mec''))'; 

timefilter{1} = {'get2dstate','($velocity<4)'};
% timefilter{2} = {'kk_getriptimes','($nripples>=1)',[],'tetfilter',tetfilter,'minthresh',5};
timefilter{2} = {'getconstimes', '($cons == 1)',[eventarea,eventtype],1,...
                   'consensus_numtets',consensus_numtets,...
                   'minthresh',minthresh,'exclusion_dur',exclusion_dur,'minvelocity',minvelocity,'maxvelocity',maxvelocity};

%----------F = createfilter('animal', animals, 'days', dayfilter,'epochs', epochfilter, 'excludetime', timefilter, 'eegtetrodes',tetfilter,'iterator', iterator);--------
F = createfilter('animal', animals, 'days', dayfilter,'epochs', epochfilter, 'excludetime', timefilter, 'eegtetrodes',tetfilter,'iterator', iterator);

%----------f = setfilteriterator(f, funcname, loadvariabl   es, options)--------
F = setfilterfunction(F, 'dfa_calcriptrigspectrogram', {'eeg', [eventarea,eventtype]},'eventtype',eventtype); 

F = runfilter(F);

%% Collect across tetrodes
    
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
    
%% Average and plot spectrograms

meanspecs = cell(1,length(tets));

for t=1:length(tets)
    meanspecs{t} = mean(F.tetout(t).S,3);
end


%plot
for t=1:length(tets)
    figure
    imagesc(F.tetout(t).t,F.tetout(t).f,meanspecs{t}',[-0.2,2.5])
    axis xy
    xlabel('Time since start of detected SWR')
    ylabel('Frequency')
    colormap('jet')
    colorbar
end
    
% save('vCA3_trigtoV_riptrigspec_run.mat','F','-v7.3')

%% plot by epochs
% 
% %plot
% for e=1:3:length(F.output{1})
%     figure
%     imagesc(F.output{1}(e).t,F.output{1}(e).f,mean(F.output{1}(e).S,3)',[-0.2,2.5])
%     axis xy
%     xlabel('Time since start of detected SWR')
%     ylabel('Frequency')
%     colormap('jet')
%     colorbar
% end
