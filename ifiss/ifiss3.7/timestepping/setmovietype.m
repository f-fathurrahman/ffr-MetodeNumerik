%SETMOVIETYPE assign movie format 
% This script is called to determine the type of video file that is
% generated and saved in the /ifiss3.7/datafiles/ directory
% Using a Windows PC, an avi video is generated and can be viewed
% using the inbuilt MATLAB function implay.m
% An avi video is also generated when running MATLAB in the cloud
% Otherwise, using Unix or Mac, mp4 files are generated

%   IFISS scriptfile: DJS; 16 November 2023.
% Copyright (c) 2023 D.J. Silvester

if ~(exist('VideoWriter')==2)
error('Oops.. movie generation is not supported (no VideoWriter function)!')
end

if (strncmp(computer,'GLNXA64',7))  %  cloud Unix architecture
avi=1; fprintf('\nGenerating .avi flow movies ... \n');    
elseif isunix
avi=0; fprintf('\nGenerating  .mp4 flow movies ... \n');
elseif ispc
avi=1; fprintf('\nGenerating  .avi flow movies ... \n');
else
error('Oops.. computer architecture not recognised!')
end
