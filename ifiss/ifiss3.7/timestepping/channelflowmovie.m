%CHANNELFLOWMOVIE channel flow solution animation
%   IFISS scriptfile: DJS; 17 November 2023.
% Copyright (c) 2023 D.J. Silvester

% fix movie format
setmovietype, help setmovietype
if avi==1,
uvid=VideoWriter('VelocityMovie.avi','Uncompressed AVI');
upvid=VideoWriter('VelocityProfileMovie.avi','Uncompressed AVI');
else
uvid=VideoWriter('VelocityMovie.mp4','MPEG-4');
uvid.FrameRate=fps; uvid.Quality=100;
upvid=VideoWriter('VelocityProfileMovie.mp4','MPEG-4');
upvid.FrameRate=fps; upvid.Quality=100;
end
% generate movie file
open(uvid);
open(upvid)
load('MovieSegments.mat');                      %  number of segments in krestart
load('VelocityScaling.mat');                    %  maximum x and y velocities
for i=1:krestart
   fname=sprintf('%s%i','MovieData',i);         %  filename of a segment
   load(fname);
   sseg=length(soltime);                        %  number of time steps in a segment
   fprintf('   Segment %i: Time steps %i to %i\n',i,(i-1)*200+1,(i-1)*200+sseg);
   channel_unsteadyflowref(qmethod,U,soltime,By,Bx,A,xy,x,y,...
                           bound,bnde,xref,i,krestart,vfreq,uvid,upvid,uxmax,uxmin);
   % compute the final divergence error
   if qmethod>1 error_div = q2div(xy,mv,[U(:,end);gzero]); end
   % compute mean quantities
   [ke,acc,meanv] = energymeanvorticity(qmethod,mv,U,Udot,soltime,By,Bx,G,xy,1,101,(i-1)*200);
   clear U Udot soltime
end
close(uvid);     %  close velocity video file
close(upvid);    %  close velocity profile video file
fprintf('All done\n')
