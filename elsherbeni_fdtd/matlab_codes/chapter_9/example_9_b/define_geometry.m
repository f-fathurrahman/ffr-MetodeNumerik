disp('defining the problem geometry');

bricks  = [];
spheres = [];

% define dielectric
bricks(1).min_x = 0;
bricks(1).min_y = 0;
bricks(1).min_z = 0;
bricks(1).max_x = 14.3e-3;
bricks(1).max_y = 25.4e-3;
bricks(1).max_z = 26.1e-3;
bricks(1).material_type = 4; 

bricks(2).min_x = -8e-3;
bricks(2).min_y = -6e-3;
bricks(2).min_z = 0;
bricks(2).max_x = 22e-3;
bricks(2).max_y = 31e-3;
bricks(2).max_z = 0;
bricks(2).material_type = 2; 

bricks(3).min_x = 0;
bricks(3).min_y = 12.2e-3;
bricks(3).min_z = 1e-3;
bricks(3).max_x = 0;
bricks(3).max_y = 13.2e-3;
bricks(3).max_z = 10e-3;
bricks(3).material_type = 2; 

