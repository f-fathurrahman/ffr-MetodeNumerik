disp('defining the problem geometry');

bricks  = [];
spheres = [];
thin_wires = [];

% define dielectric
bricks(1).min_x = -0.05;
bricks(1).min_y = -0.05;
bricks(1).min_z = 0;
bricks(1).max_x = 0.05;
bricks(1).max_y = 0.05;
bricks(1).max_z = 0.2;
bricks(1).material_type = 4; 
