disp('defining the problem geometry');

bricks  = [];
spheres = [];
thin_wires = [];

% define dielectric
bricks(1).min_x = -80e-3;
bricks(1).min_y = -80e-3;
bricks(1).min_z = -80e-3;
bricks(1).max_x = 80e-3;
bricks(1).max_y = 80e-3;
bricks(1).max_z = 80e-3;
bricks(1).material_type = 4; 
