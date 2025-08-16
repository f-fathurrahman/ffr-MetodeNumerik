clear all; close all;
%   SYNTAX 
%   examples_plyread_plywrite
%   DESCRIPTION 
%   The Polygon file format (PLY), also know as the Stanford Triangle 
%   Format, is an extremely popular and standard file format supported by
%   major CAD software packages. This script provides example usage of the 
%   files plyread.m and plywrite.m MATLAB functions featured in chapter 2 
%   of the text and written by Pascal Getreuer. These functions are
%   very important when translating between CAD and simulation software
%   platforms. 
%
%   Low-Frequency Electromagnetic Modeling for Electrical and Biological
%   Systems Using MATLAB, Sergey N. Makarov, Gregory M. Noetscher, and Ara
%   Nazarian, Wiley, New York, 2015, 1st ed.

%% Example 1: Import several PLY files using plyread.m
clear all; close all; clc;

% Select a file to import
[FileName,PathName] = uigetfile('*.ply','Select the desired PLY file to open');
file1 = cat(2,PathName,FileName);

% Invoke the plyread.m function to obtain the connectivity (tri) and vertex
% (P) matrices.
[tri1,P1] = plyread(file1,'tri');

% Visualize surface triangular elements
h1 = patch('Vertices',P1,'Faces', tri1,'facecolor','w');
hold on;
axis equal; axis tight;
view(-45,20)

% Select a second file to import
[FileName,PathName] = uigetfile('*.ply','Select the desired PLY file to open');
file2 = cat(2,PathName,FileName);

% Invoke the plyread.m function to obtain the connectivity (tri) and vertex
% (P) matrices.
[tri2,P2] = plyread(file2,'tri');

% Visualize surface triangular elements
h3 = patch('Vertices',P2,'Faces', tri2,'facecolor','b','facealpha',0.5);
hold on;
axis equal; axis tight;
view(-45,20)

%% Example 2: Write to memory an PLY file using plywrite.m
clear all; close all; clc;

% Import a file to start
[FileName,PathName] = uigetfile('*.ply','Select the desired PLY file to open');
file1 = cat(2,PathName,FileName);

% Invoke the plyread.m function to obtain the connectivity (tri) and vertex
% (P) matrices.
[tri1,P1] = plyread(file1,'tri');

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

% Write the new mesh to memory using plywrite

% Construct a structure containing face and vertex data
Data.vertex.x = P2(:,1);
Data.vertex.y = P2(:,2);
Data.vertex.z = P2(:,3);
Data.face.vertex_indices = tri1;

% Create a file to save the data
[FileName,PathName] = uiputfile('*.ply','Save PLY file as');
filesave = cat(2,PathName,FileName);

plywrite(Data,filesave);
