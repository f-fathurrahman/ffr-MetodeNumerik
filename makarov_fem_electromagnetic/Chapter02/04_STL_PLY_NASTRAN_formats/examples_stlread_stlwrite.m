%   SYNTAX 
%   examples_stlread_stlwrite
%   DESCRIPTION 
%   The STereoLithography (STL) file format is an extremely popular and
%   standard file format supported by major CAD software packages. This
%   script provides example usage of the files stlread.m and stlwrite.m
%   MATLAB functions featured in chapter 2 of the text and written by Eric
%   C. Johnson and Sven Holcombe, respectively. These functions are
%   very important when translating between CAD and simulation software
%   platforms.
%
%   Low-Frequency Electromagnetic Modeling for Electrical and Biological
%   Systems Using MATLAB, Sergey N. Makarov, Gregory M. Noetscher, and Ara
%   Nazarian, Wiley, New York, 2015, 1st ed.

%% Examples reading and writing STL files
 

%% Example 1: Import several STL files using stlread.m

% Select a file to import
[FileName,PathName] = uigetfile('*.stl','Select the desired STL file to open');
file1 = cat(2,PathName,FileName);

% Invoke the stlread.m function to obtain the connectivity (tri), vertex
% (P) and triangular element normal (Norm) matrices.
[tri1,P1,Norm1] = stlread(file1);

% Visualize surface triangular elements
h1 = patch('Vertices',P1,'Faces', tri1,'facecolor','w');
hold on;
axis equal; axis tight;
view(-45,20)

% Visualize normals
Center = (P1(tri1(:, 1),:)+P1(tri1(:, 2),:)+P1(tri1(:, 3),:))/3;
h2 = quiver3(Center(:,1),Center(:,2),Center(:,3), ...
          Norm1(:,1),Norm1(:,2),Norm1(:,3), 5, 'color', 'r');

pause

set(h2,'visible','off')

% Select a second file to import
[FileName,PathName] = uigetfile('*.stl','Select the desired STL file to open');
file2 = cat(2,PathName,FileName);

% Invoke the stlread.m function to obtain the connectivity (tri), vertex
% (P) and triangular element normal (Norm) matrices.
[tri2,P2,Norm2] = stlread(file2);

% Visualize surface triangular elements
h3 = patch('Vertices',P2,'Faces', tri2,'facecolor','b','facealpha',0.5);
hold on;
axis equal; axis tight;
view(-45,20)

%% Example 2: Write to memory an STL file using stlwrite.m
clear all; close all; clc;

% Import a file to start
[FileName,PathName] = uigetfile('*.stl','Select the desired STL file to open');
file1 = cat(2,PathName,FileName);

% Invoke the stlread.m function to obtain the connectivity (tri), vertex
% (P) and triangular element normal (Norm) matrices.
[tri1,P1,Norm1] = stlread(file1);

% Shift the nodes in the original file by a fixed distance in one direction
P2 = P1;
P2(:,1) = P1(:,1) + 600;

% Visualize surface triangular elements
h1 = patch('Vertices',P1,'Faces', tri1,'facecolor','w');
hold on;
axis equal; axis tight;
view(-45,20)

h2 = patch('Vertices',P2,'Faces', tri1,'facecolor','b');
hold on;
axis equal; axis tight;
view(-45,20)

% Write the new mesh to memory using stlwrite

% Construct a structure containing face and vertex data
FV.faces = tri1;
FV.vertices = P2;

% Create a file to save the data
[FileName,PathName] = uiputfile('*.stl','Save STL file as');
filesave = cat(2,PathName,FileName);

stlwrite(filesave,FV);
