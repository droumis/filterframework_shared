function M = makempgframe(datadir,animprefix,daynum,mpegfile,second)

% This function saves an mpeg frame from a single day's mpeg as .mat file, for use in ms_createtaskstruct.
% Particularly useful for the multiple W track where an animal may not
% visit every arm in a day. 
% INPUTS:      datadir = directory to save the mpgframe file
%               animprefix = animal's lowercase, three-letter prefix
%               mpegfile = name of the mpeg file in AnimalData/animalXX day directory
%               second = the number of seconds into the video
% OUTPUT:       M, a struct containing a MxNx3 image matrix, cdata.

% ***must start in animal's day directory***

M = mpgread(mpegfile,second*30,'truecolor');

% plot to check
figure
h = image(M.cdata);


%save mpgframe file
cd(datadir)
filename = sprintf('%smpgframe%02d.mat',animprefix,daynum);
save(filename,'M','-v7.3');
end
