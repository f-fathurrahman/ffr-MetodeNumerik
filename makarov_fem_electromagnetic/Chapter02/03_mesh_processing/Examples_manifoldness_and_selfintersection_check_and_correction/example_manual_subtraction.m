clear all; close all; warning off
%   SYNTAX 
%   example_manual_subtraction
%   DESCRIPTION 
%   This script employs manual triangle subtraction from the mesh. It is
%   useful to clean up dirty meshes. The script uses MATLAB script
%   select3D.m by Joe Conti to select triangles by right mouse click. The
%   selected triangle should be highlighted. Then, proceed to the next
%   triangle. Hit outside the mesh to remove the selected triangles and
%   stop the process. The resulting mesh will be saved in ***_mod.mat file
%
%   Low-Frequency Electromagnetic Modeling for Electrical and Biological
%   Systems Using MATLAB, Sergey N. Makarov, Gregory M. Noetscher, and Ara
%   Nazarian, Wiley, New York, 2015, 1st ed.

[FileName, PathName] = uigetfile('*.mat','Select the file');
load(FileName);
[str.MasterMeshQ str.MasterMeshIndex] = min(simpqual(P, t));
str

%   Find border triangles and border edges (if any)
[NonManifoldAttached, edges, edgesb, edgesnm] = meshedges(P, t);    

%   Plot all patches simultaneously
Color = repmat([1 0.75 0.65], size(t, 1), 1);
Color(NonManifoldAttached, :)              = repmat([0 1 0], length(NonManifoldAttached), 1);
Patch = patch('Faces', t, 'Vertices', P, 'EdgeColor', 'k', 'FaceAlpha', 1.0);
set(Patch, 'FaceColor', 'flat', 'FaceVertexCData', Color);

%   Display boundary triangles if any (only edges)
if ~isempty(edgesb)
    for m = 1:size(edgesb, 1)
        marker1(m) = line('xdata', P(edgesb(m, :),1) ,'ydata', P(edgesb(m, :),2) ,'zdata', P(edgesb(m, :),3),...
        'color', 'b', 'linewidth', 2);           
    end
end
if ~isempty(edgesnm)
    for m = 1:size(edgesnm, 1)
        marker2(m) = line('xdata', P(edgesnm(m, :),1) ,'ydata', P(edgesnm(m, :),2) ,'zdata', P(edgesnm(m, :),3),...
        'color', 'b', 'linewidth', 2);           
    end
end
axis 'equal';    axis 'tight', set(gca, 'YDir','normal');
xlabel('x, m'); ylabel('y, m'); view(-9, 19); grid on;
pause;

%% Selection loop
selection = [];
while 1   
    % p - clicked point; % v - nearest vertex; vi - index into the nearest
    % vertex; facei - index into the clicked face
    [p v vi face facei] = select3d;    
    if ~isempty(facei)
        selection = [selection facei];
        marker1 = line('xdata',p(1),'ydata',p(2),'zdata',p(3),'marker','o',...
                'erasemode','xor','markerfacecolor','k');
        marker2 = line('xdata',v(1),'ydata',v(2),'zdata',v(3),'marker','o',...
                'erasemode','xor','markerfacecolor','k');
        marker2 = line('erasemode','xor','xdata',[face(1,:) face(1,1)],...
                                         'ydata',[face(2,:) face(2,1)],...
                                         'zdata',[face(3,:) face(3,1)],'linewidth',5);
    else
        break;
    end
    pause;
end
selection
clf; t(selection, :) = [];
patch('Faces', t, 'Vertices', P, 'FaceColor', [1 0.75 0.65], 'EdgeColor', 'k', 'FaceAlpha', 1.0);
axis 'equal';    axis 'tight', set(gca, 'YDir','normal');
xlabel('x, m'); ylabel('y, m'); view(-9, 19); grid off;

[P t] = fixmesh(P, t);
NewName =  strcat(FileName(1:end-4), '_mod', '.mat');
save(NewName, 'P', 't');
