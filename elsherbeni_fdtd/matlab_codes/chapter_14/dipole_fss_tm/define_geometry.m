disp('defining the problem geometry'); 

% a brick
bricks(1).min_x =0;
bricks(1).min_y = 0;
bricks(1).min_z = -6e-3;
bricks(1).max_x = 15e-3;
bricks(1).max_y = 15e-3;
bricks(1).max_z = 0;
bricks(1).material_type = 4;

% a brick
bricks(2).min_x = 6e-3;
bricks(2).min_y = 1.5e-3;
bricks(2).min_z = 0;
bricks(2).max_x = 9e-3;
bricks(2).max_y = 13.5e-3;
bricks(2).max_z = 0.0001;
bricks(2).material_type = 2;


