function ms_lineardayprocess(directoryname,fileprefix,days, varargin)
%LINEARDAYPROCESS(directoryname,fileprefix,days, options)
%
%Runs linearizeposition for all run epochs in each day and saves the data in
%'linpos' in the directoryname folder.  See LINEARIZEPOSITION for the definitions 
%of the options.
%
%directoryname - example 'data99/user/animaldatafolder/', a folder 
%                containing processed matlab data for the animal
%
%fileprefix -    animal specific prefix for each datafile
%
%days -          a vector of experiment day numbers 
%
%                   'manual_segment' - [day1 epoch1 startindex1 endindex1 SEGMENTnumber1 ;  day2 epoch2 ...]
                        %   each row corresponds to a set of indices that
                        %   need to be manually changed -- this is after
                        %   you've verified separately that these indices
                        %   are currently being misclassified
                        %   (dfskk_thetaprecess)
                        % added kk 3.24.14
%               maxv % specify
% 
                        
welldist = 8;
maxv = 500;
manual_segment=[];

lowercasethree = '';

%set variable options
for option = 1:2:length(varargin)-1
    switch varargin{option}
        case 'lowercasethree'
            lowercasethree = varargin{option+1};
        case 'manual_segment'
            manual_segment = varargin{option+1};
        case 'maxv'
            maxv = varargin{option+1};  
        case 'welldist'
            welldist = varargin{option+1};
    end
end


days = days(:)';

for day = days
   linpos = [];
   
   dsz = '';
   if (day < 10)
      dsz = '0';
   end

   eval(['load ',directoryname,fileprefix,'task',dsz,num2str(day), '.mat']);
   eval(['task = ',lowercasethree,'task;'])
   for i = 1:length(task{day})   
      if (~isempty(task{day}{i})) && ~isfield(task{day}{i},'type')
         disp(sprintf('this epoch (d%d e%d) has no .type field!',day,i))
      elseif ((~isempty(task{day}{i})) && (strcmp(task{day}{i}.type,'run')) && (~strcmp(task{day}{i}.environment,'OpenFieldA')) && (~strcmp(task{day}{i}.environment,'OpenFieldB')))
            disp(['Day ',num2str(day), ', Epoch ',num2str(i)])
            index = [day i];
            [linpos{day}{i}.statematrix,linpos{day}{i}.segmenttable, linpos{day}{i}.trajwells, ...
                linpos{day}{i}.wellSegmentInfo, linpos{day}{i}.segmentInfo]...
                = ms_linearizeposition(directoryname,fileprefix, index, 'welldist',welldist,'maxvelocity',maxv);
            linpos{day}{i}.statematrixfields = ['time startwell endwell segment segheaddir velocity lineardist'];
            linpos{day}{i}.segmenttablefields = ['trajnum segnum segmentID'];

            % Get and save behavestate file data. This is IMP. Need it run with velocity criterion also.
            % linearizeposition only does gettraj
            % Run getbehavestate and update statematrix.traj(state) and statematrix.lindist          
%             includeStates = 6; %Means 1:6
%             [statematrix.traj, statematrix.lindist] = sj_getbehavestate(linpos,day,i, includeStates, 'minlinvelocity', 3);           
%             behavestate{day}{i}.state = linpos{day}{i}.statematrix.traj;
%             behavestate{day}{i}.lindist = linpos{day}{i}.statematrix.lindist;
          
      end
   end
   
   eval([lowercasethree,'linpos = linpos;']);
   eval(['save ',directoryname,fileprefix,'linpos',dsz,num2str(day),' ',lowercasethree,'linpos']);
   
   % Save behavestatefile
%    statefile = sprintf('%s/%sbehavestate%02d.mat', directoryname, fileprefix, day);
%    save(statefile,'behavestate');    
end

