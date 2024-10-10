clear all; close all;

walldatain = fopen('walldata.dat');
% This will read all lines
lines_txt1 = textscan(walldatain, '%s', 'delimiter', '\n', 'CommentStyle', '%', 'MultipleDelimsAsOne',1);
fclose(walldatain);
lines_txt = lines_txt1{1};

iline = 1;
wallNum = str2double(lines_txt{iline});

for i=1:wallNum
  %
  iline = iline + 1;
  numcorn = str2double(lines_txt{iline});
  wallSide(i) = numcorn;
  %
  iline = iline + 1;
  isfree = str2double(lines_txt{iline});
  wallisfree(i,1) = isfree;
  %
  for j = 1:wallSide(i)
    iline = iline + 1;
    wallpos = str2num(lines_txt{iline});
    wallX(j,i) = wallpos(1);
    wallY(j,i) = wallpos(2);
  end
end

