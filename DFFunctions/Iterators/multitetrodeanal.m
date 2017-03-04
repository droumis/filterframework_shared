function f = multitetrodeanal(f)
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
    totalepochs = [];
    for g = 1:length(f(an).epochs)
        totalepochs = [totalepochs; f(an).epochs{g}];
    end
    totaldays = unique(totalepochs(:,1)); %get all of the days across groups
    
    %load all the variables that the function requires except the eeg
    loadstring = [];
    for i = 1:length(f(an).function.loadvariables)
        
            eval([f(an).function.loadvariables{i},' = loaddatastruct(animaldir, animalprefix, f(an).function.loadvariables{i}, totaldays);']);
            loadstring = [loadstring, f(an).function.loadvariables{i},','];
    end
    foptions = f(an).function.options;
    
    %iterate through the epochs within each data group
    for g = 1:length(f(an).epochs)
        
        for e = 1:size(f(an).epochs{g},1)
            
            % load the eeg data for all the specified ntrodes
            for i = 1:length(f(an).function.loadvariables)
               if (iseegvar(f(an).function.loadvariables{i}))
                    loadday=f(an).epochs{g}(e,1);
                    loadepoch=f(an).epochs{g}(e,2);
                    loadntrodes = f(an).eegdata{g}{e};
                    eval([f(an).function.loadvariables{i},' = loadeegstruct(animaldir, animalprefix, f(an).function.loadvariables{i}, loadday, loadepoch, loadntrodes);']);
               end
            end                 
            
            
            
            indices = f(an).eegdata{g}{e};
            if size(indices,2) > 1
                error(['Data must only include tetrodes'])
            end
            numindices = size(indices,1);
            indices = [repmat(f(an).epochs{g}(e,:),[numindices 1]) indices];
            excludeperiods = f(an).excludetime{g}{e};
            % run the specified filter function on this set of animal/epoch/ntrodes
            eval(['fout = ',f(an).function.name,'(indices,excludeperiods,', loadstring, 'foptions{:});']);

            %save the function output in the filter variable.  Allows numeric or struct outputs
            if isstruct(fout)
                f(an).output{g}(e) = fout;
            elseif isnumeric(fout)          
                if (length(f(an).output) < g)
                    f(an).output{g} = [];
                end
                f(an).output{g} = stack(f(an).output{g}, fout);             
            else
                error(['In calling ', f(an).function.name, ': Function output must be either numeric or a structure']);
            end
            
        end
    end
end