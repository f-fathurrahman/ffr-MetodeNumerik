disp('defining the problem geometry'); 

% ground p1ane
bricks(1).min_x = 0;
bricks(1).min_y = 0;
bricks(1).min_z = 0;
bricks(1).max_x = 60e-3;
bricks(1).max_y = 40e-3;
bricks(1).max_z = 0;
bricks(1).material_type = 2;

% substrate
bricks(2).min_x = 0;
bricks(2).min_y = 0;
bricks(2).min_z = 0;
bricks(2).max_x = 60e-3;
bricks(2).max_y = 40e-3;
bricks(2).max_z = 1.9e-3;
bricks(2).material_type = 4; % er = 2.2

% patch
bricks(3).min_x = 2e-3;
bricks(3).min_y = 10e-3;
bricks(3).min_z = 1.9e-3;
bricks(3).max_x = 58e-3;
bricks(3).max_y = 30e-3;
bricks(3).max_z = 1.9e-3;
bricks(3).material_type = 2;


