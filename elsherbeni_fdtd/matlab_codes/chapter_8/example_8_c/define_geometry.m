disp('defining the problem geometry');

bricks  = [];
spheres = [];

% define a substrate
bricks(1).min_x = 0;
bricks(1).min_y = 0;
bricks(1).min_z = 0;
bricks(1).max_x = 60*dx;
bricks(1).max_y = 60*dy;
bricks(1).max_z = 6*dz;
bricks(1).material_type = 4; 

% define a PEC plate 
bricks(2).min_x = 24*dx;
bricks(2).min_y = 0;
bricks(2).min_z = 6*dz;
bricks(2).max_x = 36*dx;
bricks(2).max_y = 60*dy;
bricks(2).max_z = 6.02*dz;
bricks(2).material_type = 2; 
