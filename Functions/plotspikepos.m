function h = ms_plotspikepos(spikes, pos, ind, color, varargin)
%Modification of function plotspikepos(spikes, pos, ind)
%  - plots positions in grey and spikes as the specified color, which must be
%  a vector of 3 numbers ([r g b])
%  - looks up spike times in pos data to get spike position
%  - have to load spikes and pos into workspace eithe manually or with a
%  calling script
%  - ind specifies [day epoch tetrode cell]

day = ind(1);
ep = ind(2);
tet = ind(3);
cell = ind(4);
p = pos{day}{ep}.data;
s = spikes{day}{ep}{tet}{cell}.data;
halftocalc = 0; %default 0 if plotting whole epoch
suppresspos = 0; %default don't suppress position plotting
spikemarker = '.';
jitter = 0; 

for option = 1:2:length(varargin)-1
    switch varargin{option}
        case 'half'
            halftocalc = varargin{option+1};
            case 'suppresspos'
            suppresspos = varargin{option+1};
        case 'spikemarker'
            spikemarker = varargin{option+1};
        case 'jitter'
            jitter = varargin{option+1};
        otherwise
            error(['Option ''', varargin{option}, ''' not defined']);
    end
end


% find position at all spike times
X = p(lookup(s(:,1),p(:,1)),2);
Y = p(lookup(s(:,1),p(:,1)),3);


if halftocalc>0
midtime=pos{day}{ep}.data(1,1)+((pos{day}{ep}.data(end,1)-pos{day}{ep}.data(1,1))/2);
midposindex = lookup(midtime,pos{day}{ep}.data(:,1));
midspikeindex = lookup(midtime,spikes{day}{ep}{tet}{cell}.data(:,1));
end

if halftocalc==1 %Fix with new X and Y
%     figure
    if suppresspos==0
    h1 = plot(p(1:midposindex,2), p(1:midposindex,3), '.');
    set(h1, 'color', [0.8 0.8 0.8]);
    set(h1, 'MarkerSize', 6);
    end
    
    hold on
    h2 = plot(s(1:midspikeindex,2), s(1:midspikeindex,3), '.','color', color);
    set(h2, 'MarkerSize', 10);
elseif halftocalc==2
%     figure
if suppresspos==0
    h1 = plot(p(midposindex:end,2), p(midposindex:end,3), '.');
    set(h1, 'color', [0.8 0.8 0.8]);
    set(h1, 'MarkerSize', 6);
end
    
    hold on
    h2 = plot(s(midspikeindex:end,2), s(midspikeindex:end,3), '.','color', color);
    set(h2, 'MarkerSize', 10);
else
    
%     figure
if suppresspos==0
    h1 = plot(p(:,2), p(:,3), '.');
    set(h1, 'color', [0.8 0.8 0.8]);
    set(h1, 'MarkerSize', 6);
end
    if jitter~=0
         hold on
    h2 = plot(X+jitter, Y, spikemarker,'color', color);
    set(h2, 'MarkerSize', 6,'LineWidth',1);  %size 10 if dots, 4 if circles
    else
    hold on
    h2 = plot(X, Y, spikemarker,'color', color);
    set(h2, 'MarkerSize', 6,'LineWidth',1);
    end
end



