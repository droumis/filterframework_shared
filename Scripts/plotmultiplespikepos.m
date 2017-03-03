% MS 2016
% For or every cell on a given tetrode and day, makes a figure with per-epoch subplots of spike position 

%Can plot all of one day or a subset of cells from 1 animal
plotday = 1;
plotripplemodcells = 0;
savefigs = 0; %automtatically save figures?

if plotday
DAYS = 11;
else
DAYS = unique(cellindices(:,1)); % if you have sig ripple mod cell information saved somewhere
end

%% 

for day = DAYS'
animdir = '/opt/data40/mari/Eli/';
animprefix = 'eli';
%load necessary variables for function hj_plotspikepos = pos, spikes
pos=loaddatastruct(animdir,animprefix,'pos',day);
spikes=loaddatastruct(animdir,animprefix,'spikes',day);


%record number of epochs for each day as a matrix
epochs=1:length(spikes{day});
if floor(epochs(end)/2)==epochs(end)/2
subplotsize = epochs(end)/2;
else
    subplotsize = floor(epochs(end)/2) + 1;
end

%record numbers of cells for each clustered tetrode as a matrix
if plotday
    tetcells = [];
    for e = epochs
        tetlist = find(cellfun(@isempty,spikes{day}{e})==0);
        for t = tetlist
            cellIDs = find(cellfun(@isempty,spikes{day}{e}{t})==0);
            tetcells = [tetcells; repmat(t,length(cellIDs),1) cellIDs'];
        end
    end
    tetcells = unique(tetcells(:,:),'rows');
end

if plotripplemodcells
    tetcells = cellindices((cellindices(:,1)==day),2:3);
end

% Fab day 6: [25 1;26 3;27 1;27 2;27 3;27 4;28 1;28 2;28 4;22 1;1 2]
% Fab day 8: [1 2;22 1;26 3;27 1;27 2;27 3;27 4;28 1]

for cell = 1:size(tetcells,1)
    %set spike location marker color in RGB
color= [rand(1) rand(1) rand(1)];
    h = figure;
    hold on
    for d = day
        for epoch = epochs
            %find the corresponding number of cells from matrix tetcells
            %find the number of epochs for the given day
            ind = [d epoch tetcells(cell,1) tetcells(cell,2)]
            if ~isempty(spikes{ind(1)}{ind(2)}{ind(3)}{ind(4)}) && ~isempty(spikes{ind(1)}{ind(2)}{ind(3)}{ind(4)}.data)
            subplot(2,subplotsize,epoch)
            ms_plotspikepos(spikes, pos, ind, color)
            title(['D' num2str(d) 'E' num2str(epoch) 'T' num2str(tetcells(cell,1)) 'C' num2str(tetcells(cell,2))])
            end
        end
    end
    if savefigs
    filename = sprintf('%sspikeposition/allspikepos_%d-%d-%d.fig',animdir,d,tetcells(cell,1),tetcells(cell,2));
        saveas(h,filename)
        clear color
    end
end
end