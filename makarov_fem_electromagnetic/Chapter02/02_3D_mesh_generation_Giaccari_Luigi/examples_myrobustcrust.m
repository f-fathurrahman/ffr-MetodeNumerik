clear all; close all;
%   SYNTAX 
%   examples_myrobustcrust
%   DESCRIPTION 
%   This script provides examples for using the MATLAB function
%   MyRobustCrust, written by Dr. Luigi Giaccari and featured in
%   Chapter 2 of the text. 

%   Example 1 demonstrates use of the script by constructing a number of
%   random points on a unit sphere, constructing their triangulation, and
%   plotting the results. Example 2 again demonstrates use of the script
%   by importing in a pre-existing point cloud, constructing the
%   corresponding triangulation, and plotting the results.
%
%   Low-Frequency Electromagnetic Modeling for Electrical and Biological
%   Systems Using MATLAB, Sergey N. Makarov, Gregory M. Noetscher, and Ara
%   Nazarian, Wiley, New York, 2015, 1st ed.

%% Example 1: 
%   Generating a mesh of a unit sphere
clear all; close all; clc;

% Generate a number of points on a unit sphere
n=1000;                 
theta=pi*rand(1,n);
phi=2*pi*rand(1,n);
[x,y,z]=sph2cart(theta,phi,ones(1,n));

% Construct a single Points array of x-, y-, and z-coordinates
P = [x' y' z'];

% Invoke MyRobustCrust to generate the connectivity matrix
tri = MyRobustCrust(P);

%Visualize the result
h = patch('Vertices',P,'Faces', tri,'facecolor',[0.8, 0.5, 0.3]);
axis equal

%% Example 2: 
%   Generate mesh triangulation from a file containing pre-existing data
clear all; close all; clc;

% Load point cloud data
load('Points.mat')

% Invoke MyRobustCrust to generate the connectivity matrix
tri = MyRobustCrust(P);

%Visualize the results - point cloud data on the left, triangulation on the
%right

h = patch('Vertices',P,'Faces', tri,'facecolor',[0.8, 0.2, 0.3]);
subplot(1,2,1)
plot3(P(:,1), P(:,2), P(:,3), '.')
axis equal; axis tight;
view(-43,18)

subplot(1,2,2)
h = patch('Vertices',P,'Faces', tri,'facecolor',[0.8, 0.2, 0.3]);
axis equal; axis tight;
view(-43,18)