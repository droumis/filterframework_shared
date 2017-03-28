function f = singleeeganal(f, varargin)
% f = singleeeganal(f)
% Iterator for a filter object.  Calls the function designated in
% f().function.name, after loading the eeg variables designated as strings in
% f().function.loadvariables{:}.  Also the function call appends any
% options in the f().function.options{} cell array.
%
% Each function call is for one tetrode, and it is assumed that
% the function's first input is the index to the tetrode ([day epoch tetrode]).
% The second input is a list of exclusion periods [starttime endtime].
% The next inputs are the load variables.  Note that eeg load variables,
% specified in iseegvar(), are loaded individually for each tetrode.
% The final inputs are the options.
% out = fname(index, excludeperiods, var1, var2, ..., option1, option2,...).
%
% The output of the call function can either be a 1 by N vector, or a structure.
% The outputs are stored in f().output, grouped using the same groupings as
% in the filter.

%load all the variables that the function requires
foptions = f.function.options;
%find all unique epochs to analyze for the current animal
animaldir = f.animal{2};
animalprefix = f.animal{3};
totalepochs = [];
lendayeps = length(f.epochs{1});
dayeps = f.epochs{1};
days = unique(dayeps(:,1)); %get all of the days across groups
%iterate through the epochs within each data group
lendays = length(days);
for iday = 1:lendays
    day = days(iday);
    eps = dayeps(find(day==dayeps(:,1)),2); %get epochs for this daY
    lenieps = length(eps);
    for iep = 1:lenieps
        ep = eps(iep);
        lenitets = size(f.eegdata{iday}{iep},1);
        for itet = 1:lenitets
            tet = f.eegdata{iday}{iep}(itet);
            tmpindex = [day ep tet];
            excludeperiods = f.excludetime{iday}{iep};
            loadstring = [];
            lenvars = length(f.function.loadvariables);
            for ivar = 1:lenvars
                if (iseegvar(f.function.loadvariables{ivar}))
                    eval([f.function.loadvariables{ivar},' = loadeegstruct(animaldir, animalprefix, f.function.loadvariables{ivar}, tmpindex(1), tmpindex(2), tmpindex(3:end));']);
                    loadstring = [loadstring, f.function.loadvariables{ivar},','];
                end
            end
            %run the designated function: fout = fname(tmpindex, var1, var2, ..., option1, option2, ...)
            eval(['f.output{' num2str(day) '}{' num2str(ep) '}{' num2str(tet) '} = ',f.function.name,'(tmpindex,excludeperiods,', loadstring, 'foptions{:});']);
        end
    end
end
