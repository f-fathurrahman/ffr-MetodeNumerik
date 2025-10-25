%CAVITYFLOWMOVIE cavity flow solution animation
%   IFISS scriptfile: DJS; 16 November 2023.
% Copyright (c) 2023 D.J. Silvester

% fix movie format
setmovietype, help setmovietype
if avi==1,
uvid=VideoWriter('VelocityMovie.avi','Uncompressed AVI');
vvid=VideoWriter('VorticityMovie.avi','Uncompressed AVI');
uxpvid=VideoWriter('xVelocityProfileMovie.avi','Uncompressed AVI');
uypvid=VideoWriter('yVelocityProfileMovie.avi','Uncompressed AVI');
else
uvid=VideoWriter('VelocityMovie.mp4','MPEG-4');
uvid.FrameRate=fps; uvid.Quality=100;
vvid=VideoWriter('VorticityMovie.mp4','MPEG-4');
vvid.FrameRate=fps; vvid.Quality=100;
uxpvid=VideoWriter('xVelocityProfile.mp4','MPEG-4');
uxpvid.FrameRate=fps; uxpvid.Quality=100;
uypvid=VideoWriter('yVelocityProfile.mp4','MPEG-4');
uypvid.FrameRate=fps; uypvid.Quality=100;
end
% generate movie file
open(uvid);
open(vvid);
open(uxpvid);
open(uypvid);
load('MovieSegments.mat');                 %  number of segments in krestart
load('VelocityScaling.mat');               %  maximum x and y velocities
for i=1:krestart
   fname=sprintf('%s%i','MovieData',i);    %  filename of a segment
   load(fname);
   sseg=length(soltime);                   %  number of time steps in a segment
   fprintf('   Segment %i: Time steps %i to %i\n',i,(i-1)*200+1,(i-1)*200+sseg);
   square_unsteadyflowref(qmethod,U,soltime,mv,By,Bx,A,G,xy,xyp,x,y,bound,...
    xref,yref,i,krestart,vfreq,uvid,vvid,uxpvid,uypvid,uxmax,uxmin,uymax,uymin)
   % compute the final divergence error
   if qmethod>1 error_div = q2div(xy,mv,[U(:,end);gzero]); end
   % compute mean quantities
   [ke,acc,meanv] = energymeanvorticity(qmethod,mv,U,Udot,soltime,By,Bx,G,xy,1,101,(i-1)*200);
   clear U Udot soltime
end
close(uvid);     %  close velocity video file
close(vvid);     %  close vorticity video file
close(uxpvid);   %  close x-velocity profile video file
close(uypvid);   %  close y-velocity profile video file
fprintf('All done\n')
