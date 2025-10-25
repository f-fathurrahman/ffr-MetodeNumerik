%STEPBOUSSMOVIE Boussinesq step flow solution animation
%   IFISS scriptfile: DJS; 16 November 2023.
% Copyright (c) 2023 D.J. Silvester

% fix movie format
setmovietype, help setmovietype
if avi==1,
uvid=VideoWriter('VelocityMovie.avi','Uncompressed AVI');
tvid=VideoWriter('TemperatureMovie.avi','Uncompressed AVI');
vvid=VideoWriter('VorticityMovie.avi','Uncompressed AVI');
upvid=VideoWriter('VelocityProfileMovie.avi','Uncompressed AVI');
else
uvid=VideoWriter('VelocityMovie.mp4','MPEG-4');
uvid.FrameRate=fps; uvid.Quality=100;
tvid=VideoWriter('TemperatureMovie.mp4','MPEG-4');
tvid.FrameRate=fps; tvid.Quality=100;
vvid=VideoWriter('VorticityMovie.mp4','MPEG-4');
vvid.FrameRate=fps; vvid.Quality=100;
upvid=VideoWriter('VelocityProfileMovie.mp4','MPEG-4');
upvid.FrameRate=fps; upvid.Quality=100;
end

open(uvid);                                %  open velocity video file
open(tvid);                                %  open temperature video file
open(vvid);                                %  open vorticity video file
open(upvid);                               %  open velocity profile video file
load('MovieSegments.mat');                 %  number of segments in krestart
load('VelocityScaling.mat');               %  maximum x and y velocities
tout=0; unpack_boussdata
for i=1:krestart
   fr=[];                                 %  flow rates through the intersection
   tfr=[];                                %  discrete times of flow rates
   fname=sprintf('%s%i','MovieData',i);   %  filename of a segment
   load(fname);
   sseg=length(soltime);
   fprintf('   Segment %i: Time steps %i to %i\n',i,(i-1)*200+1,(i-1)*200+sseg);
      step_unsteadytempref(T,soltime,xyt,grid(1).x,grid(1).y,vfreq,i,krestart,tvid);
      step_unsteadyflowref(2,mv2,U,soltime,Av,BBy,BBx,Qv,xyv,xyp,grid(1).x, ...
            grid(1).y,bnd_d,i,krestart,vfreq,xref,uvid,vvid,upvid,uxmax,0);
 %  figure(2000);
 %  plot(tfr,fr,'.-b'); hold on
 %  xlabel('time'); ylabel('flow rate');
 %  title('Flow rate through the middle of the hot segment');
   error_div = q2div(xyv,grid(1).mv2,[U(:,end);gzero]);
   clear U T soltime
end
close(uvid);     %  close velocity video file
close(vvid);     %  close vorticity video file
close(tvid);     %  close temperature video file
close(upvid);    %  close velocity profile video file
fprintf('All done\n')
