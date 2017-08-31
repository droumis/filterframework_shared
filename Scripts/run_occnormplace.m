%% Calculates and plots occupancy normalized place field using function occnormplace

% Load data and set up parameters
animprefix = 'fab';
detc = [16 3 27 3];
day = detc(1);
ep = detc(2);
tet = detc(3);
cell = detc(4);
% plot for only half the epoch?
half = 0; % 1 = 1st half, 2 = 2nd half, 0 = plot whole epoch

pos=loaddatastruct('/opt/data40/mari/Fab/','fab','pos',detc(1));
spikes=loaddatastruct('/opt/data40/mari/Fab/','fab','spikes',detc(1));
% pos=loaddatastruct('/opt/data40/mari/Eli/','eli','pos',detc(1));
% spikes=loaddatastruct('/opt/data40/mari/Eli/','eli','spikes',detc(1));

% load(sprintf('%spos%d.mat',animprefix,day));
% load(sprintf('%sspikes%d.mat',animprefix,day));

midtime=pos{day}{ep}.data(1,1)+((pos{day}{ep}.data(end,1)-pos{day}{ep}.data(1,1))/2);
midposindex = lookup(midtime,pos{day}{ep}.data(:,1));
midspikeindex = lookup(midtime,spikes{day}{ep}{tet}{cell}.data(:,1));

allposdata = pos{day}{ep}.data;
if half == 1
    xypos = allposdata(1:midposindex,2:3);
    allxypos = allposdata(:,2:3);
    spikepos = spikes{day}{ep}{tet}{cell}.data(1:midspikeindex,2:3);
elseif half == 2
    xypos = allposdata(midposindex:end,2:3);
    allxypos = allposdata(:,2:3);
    spikepos = spikes{day}{ep}{tet}{cell}.data(midspikeindex:end,2:3);
else
    xypos = allposdata(:,2:3);
    spikepos = spikes{day}{ep}{tet}{cell}.data(:,2:3);
end
timestep = allposdata(2,1) - allposdata(1,1);
std = 4; %std of smoothing Gaussian, in cm
binsize = 1; %cm

out = occnormplace(xypos,spikepos,binsize,timestep,std)

%% Plot occnormplace output

output = out;
occnormflag =0;
cutoff_percentile = 5;

%%% PLOTTING PARAMETERS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
data_toplot = output.eventrate_smoothed;

% if occnormflag
%     data_toplot = output.eventrate_smoothed;
% else
%     data_toplot = output.events_smoothed;
% end

%%% PLOTTING PARAMETERS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure
hold on
% Keeps the proportions of the data
% axis image;
% Labels the sides as 'cm'
xlabel('cm');
ylabel('cm');
% The tick marks are  apart
maxx = max(output.binx);
maxy = max(output.biny);
minx = min(output.binx);
miny = min(output.biny);
set(gca,'xlim',[0 maxx+5]);
% set(gca, 'XTick', [0:5:(maxx-minx)+5]);
set(gca,'ylim',[0 maxy+5]);
% set(gca, 'YTick', [0:5:(maxy-miny)+5]);
% colormap stuff
nc = 1024;
% if occnormflag
    cmap = jet(nc);
% else
%     cmap = hot(nc);
% end

cmap(1,:) = 1;      % lowest value gets a white background
colormap(cmap);

% mask the areas in which the animal actually travelled
if 1
    occ_toplot = output.occupancy_smoothed;
    cutoffvalue = prctile(10*output.occupancy_smoothed(:),cutoff_percentile);    
    MASK = (output.occupancy_smoothed > cutoffvalue);  
    data_toplot(~MASK) = -0.02;
%     occ_toplot(~MASK) = -0.02;
%     occ_toplot(MASK) = 0;
    % if it's the smoothed eventrate, then take MASK areas that have 0
    % values and make them a small positive value so they show up in the
    % plot not as white
    if occnormflag
        whitevals = (data_toplot <= 0.02);
        data_toplot(whitevals & MASK) = 0.01;
    end
end
% 
maxvalue = max(data_toplot(:))
% maxvalue = 7;

% PLOT %%%%%%%%%%%%%%%%%%%%%%
imagesc(output.biny,output.binx,data_toplot') %'CDataMapping','scaled');
% imagesc(output.epochbiny,output.epochbinx,occ_toplot')
b = colorbar('EastOutside');
caxis([-0.01 maxvalue]); %maxvalue
p = get(b, 'Position');

freezeColors

set(gca,'YDir','normal')
set(gca,'XDir','normal')

axis tight

