function f = multitetrodeanal(f, varargin)
% f = multicellanal(f)
% Iterator for a filter object.  Calls the function designated in
% f().function.name, after loading the variables designated as strings in
% f().function.loadvariables{:}.  Also the function call appends any
% options in the f().function.options{} cell array.
%
% Each function call is for one epoch, and it is assumed that
% the function's first input is a list of indices to the cell or tetrode ([day epoch tetrode
% cell]).  The second input is a list of exclusion periods [starttime endtime].
% The next inputs are the load variables, and the final inputs are the options.
% out = fname(index, excludeperiods, var1, var2, ..., option1, option2,...).
%
% The output of the call function can either be numeric, or a structure.
% If the output if numeric, data is appended along the first dimenion for
% the group.
% The outputs are stored in f().output, grouped using the same groupings as
% in the filter.

%iterate through all animals
for an = 1:length(f)
    %find all unique epochs to analyze for the current animal
    animaldir = f(an).animal{2};
    animalprefix = f(an).animal{3};
    totaldayepochs = [];
    for iepochs = 1:length(f(an).epochs)
        totaldayepochs = [totaldayepochs; f(an).epochs{1}];
    end
    totaldays = unique(f(an).epochs{1}(:,1)); %get all of the days across groups
    totalepochs = totaldayepochs(:,2);
    
    %load all the variables that the function requires except the eeg
    loadstring = [];
    for i = 1:length(f(an).function.loadvariables)
%         disp(f(an).function.loadvariables{i})
        eval([f(an).function.loadvariables{i},' = loaddatastruct(animaldir, animalprefix, f(an).function.loadvariables{i}, totaldays);']);
        loadstring = [loadstring, f(an).function.loadvariables{i},','];
    end
    foptions = f(an).function.options;
    
    %iterate through the epochs within each data group
    for idayep = 1:length(f(an).epochs{1}(:,1))
        day = f(an).epochs{1}(idayep,1);
        epoch = f(an).epochs{1}(idayep,2);
        % load the eeg data for all the specified ntrodes
        for i = 1:length(f(an).function.loadvariables)
            if (iseegvar(f(an).function.loadvariables{i}))
                ntrodes = f(an).eegdata{1}{idayep};
%                 disp(f(an).function.loadvariables{i})
                eval([f(an).function.loadvariables{i},' = loadeegstruct(animaldir, animalprefix, f(an).function.loadvariables{i}, day, epoch, ntrodes);']);
            end
        end
        ntrodes = f(an).eegdata{1}{idayep};
        if size(ntrodes,2) > 1
            error(['Data array must only include ntrodes vec'])
        end
        numntrodes = size(ntrodes,1);
        indices = [repmat(f(an).epochs{1}(idayep,:),[numntrodes 1]) ntrodes];
        excludeperiods = f(an).excludetime{1}{idayep};
        if isempty(indices)
            fprintf(sprintf('no data for %s Day%d ep%d.. skipping \n', animalprefix, day, epoch));
            continue
        end
        % run the specified filter function on this set of animal/epoch/ntrodes
        eval(['fout = ',f(an).function.name,'(indices,excludeperiods,' loadstring, 'foptions{:});']);
        if isempty(fout.data)
            continue
        end
        %save the function output in the filter variable.  Allows numeric or struct outputs
        if isstruct(fout)
            f(an).output{day}(epoch) = fout;
            
%             elseif isnumeric(fout)
%                 if (length(f(an).output) < day)
%                     f(an).output{day} = [];
%                 end
%                 f(an).output{day} = stack(f(an).output{day}, fout);
        else
            error(['In calling ', f(an).function.name, ': Function output must be either numeric or a structure']);
        end
        
    end
end
end