clear all; close all;
%   SYNTAX 
%   examples_nastran_write
%   DESCRIPTION 
%   The NASTRAN format is the oldest ASCII format for saving triangular
%   meshes. The present script converts a P, t mesh data from MATLAB to a
%   NASTRAN file.
%
%   Low-Frequency Electromagnetic Modeling for Electrical and Biological
%   Systems Using MATLAB, Sergey N. Makarov, Gregory M. Noetscher, and Ara
%   Nazarian, Wiley, New York, 2015, 1st ed.

FileName = uigetfile('*.mat','Select the mesh file to open');
S = load(FileName, '-mat');
P = S.P; t = S.t; 
fileID = fopen(strcat(FileName(1:end-4), '.nas'), 'w');

for m = 1:size(P, 1)
    fprintf(fileID, '%-8s%-8s%-8s%8s%8s%8s\n', 'GRID', num2str(m, 5), '', ...
    num2str(P(m, 1), '%6.1f'), num2str(P(m, 2), '%6.1f'), num2str(P(m, 3), '%6.1f'));
end
for m = 1:size(t, 1)
     fprintf(fileID, '%-8s%-8s%-8s%8s%8s%8s\n', 'CTRIA3', num2str(m, 5), '1', ...
     num2str(t(m, 1), '%6d'), num2str(t(m, 2), '%6d'), num2str(t(m, 3), '%6d'));   
end

fclose(fileID);

PP = P';
tt = t';
Xshape = reshape(PP(1, tt(1:3,:)),[3,size(tt,2)]);
Yshape = reshape(PP(2, tt(1:3,:)),[3,size(tt,2)]);
Zshape = reshape(PP(3, tt(1:3,:)),[3,size(tt,2)]);
patch(Xshape, Yshape, Zshape, [1 0.75 0.65], 'EdgeColor', 'k', 'FaceAlpha', 1.0,'FaceColor', [1 0.75 0.65]);
axis 'equal';  axis 'tight';
xlabel('x, mm'); ylabel('y, mm'); zlabel('z, mm');