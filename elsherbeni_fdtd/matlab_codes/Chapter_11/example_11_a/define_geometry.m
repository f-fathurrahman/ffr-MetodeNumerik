disp('defining the problem geometry');

bricks  = [];
spheres = [];
thin_wires = [];

% define a sphere 
spheres(1).radius = 100e-3;
spheres(1).center_x = 0;
spheres(1).center_y = 0;
spheres(1).center_z = 0;
spheres(1).material_type = 4;
