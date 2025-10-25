%CAVITYBOUSSMOVIE  Boussinesq cavity flow solution animation
%   IFISS scriptfile: DJS; 16 November 2023.
% Copyright (c) 2023 D.J. Silvester

% fix movie format
setmovietype, help setmovietype
if avi==1,
uvid=VideoWriter('VelocityMovie.avi','Uncompressed AVI');
tvid=VideoWriter('TemperatureMovie.avi','Uncompressed AVI');
switch(hty)                             %  cavity problem switch
case(1)   %  Rayleigh-Benard
   upvid=VideoWriter('VelocityProfileMovie.avi','Uncompressed AVI');
case(2)   %  laterally heated
   upxvid=VideoWriter('xVelocityProfileMovie.avi','Uncompressed AVI');
   upyvid=VideoWriter('yVelocityProfileMovie.avi','Uncompressed AVI');
end
else
uvid=VideoWriter('VelocityMovie.mp4','MPEG-4');
uvid.FrameRate=fps; uvid.Quality=100;
tvid=VideoWriter('TemperatureMovie.mp4','MPEG-4');
tvid.FrameRate=fps; tvid.Quality=100;
switch(hty)                             %  cavity problem switch
   case(1)   %  Rayleigh-Benard
      upvid=VideoWriter('VelocityProfileMovie.mp4','MPEG-4');
      upvid.FrameRate=fps; upvid.Quality=100;
   case(2)   %  laterally heated
      upxvid=VideoWriter('xVelocityProfileMovie.mp4','MPEG-4');
      upxvid.FrameRate=fps; upxvid.Quality=100;
      upyvid=VideoWriter('yVelocityProfileMovie.mp4','MPEG-4');
      upyvid.FrameRate=fps; upyvid.Quality=100;
      end
end
open(uvid);                       %  open velocity video file
open(tvid);                       %  open temperature video file
switch(hty)
   case(1)
      open(upvid);                %  open velocity profile video file
      upxvid=[]; upyvid=[];
   case(2)
      upvid=[];
      open(upxvid);               %  open x velocity profile video file
      open(upyvid);               %  open y velocity profile video file
end
load('MovieSegments.mat');        %  number of segments in krestart
load('VelocityScaling.mat');      %  maximum x and y velocities
tout=0; unpack_boussdata
fr=[];                            %  flow rates through the intersection
frx=[];                           %  x-flow rate
fry=[];                           %  y-flow rate
tfr=[];                           %  discrete times of flow rates
for i=1:krestart
   fname=sprintf('%s%i','MovieData',i);   %  filename of a segment
   load(fname);
   sseg=length(soltime);                  %  number of time steps in a segment
   fprintf('   Segment %i: Time steps %i to %i\n',i,(i-1)*200+1,(i-1)*200+sseg);
   [fr,frx,fry,tfr]=cavity_unsteadyflowref(2,mv2,U,T,soltime,Av,BBy,BBx,Qv,...
        xyv,grid(1).x,grid(1).y,grid(1).bnd_dn2,npx,npy,net,xref,yref,hty,i,krestart,vfreq,...
        uvid,tvid,upvid,upxvid,upyvid,uxmax,uymax,fr,frx,fry,tfr);
   % compute the divergence error
   error_div = q2div(xyv,grid(1).mv2,[U(:,end);gzero]);
   clear U T soltime
end
switch(hty)
   case(1)
      fprintf('Maximum flow rate at y=%7.4f is %6.4e\n',yref,max(fr));
   case(2)
      fprintf('Maximum flow rate at x=%7.4f is %6.4e\n',xref,max(frx));
      fprintf('Maximum flow rate at y=%7.4f is %6.4e\n',yref,max(fry));
end
close(uvid);
close(tvid);
switch(hty)
   case(1)
      close(upvid);
   case(2)
      close(upxvid);
      close(upyvid);
end
fprintf('All done\n')
