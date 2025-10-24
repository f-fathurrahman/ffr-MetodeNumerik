disp('defining the problem geometry');

bricks  = [];
spheres = [];

% define a brick with material type 4
bricks(1).min_x = 0;
bricks(1).min_y = 0;
bricks(1).min_z = 0;
bricks(1).max_x = 10e-3;
bricks(1).max_y = 18e-3;
bricks(1).max_z = 1e-3;
bricks(1).material_type = 4; 

% define a brick with material type 2
bricks(2).min_x = 4e-3;
bricks(2).min_y = 0;
bricks(2).min_z = 1e-3;
bricks(2).max_x = 5.8e-3;
bricks(2).max_y = 4e-3;
bricks(2).max_z = 1e-3;
bricks(2).material_type = 2; 

% define a brick with material type 2
bricks(3).min_x = 4.4e-3;
bricks(3).min_y = 4e-3;
bricks(3).min_z = 1e-3;
bricks(3).max_x = 5.4e-3;
bricks(3).max_y = 14e-3;
bricks(3).max_z = 1e-3;
bricks(3).material_type = 2; 

% define a brick with material type 2
bricks(4).min_x = 4.8e-3;
bricks(4).min_y = 14e-3;
bricks(4).min_z = 1e-3;
bricks(4).max_x = 5.2e-3;
bricks(4).max_y = 18e-3;
bricks(4).max_z = 1e-3;
bricks(4).material_type = 2; 

% define absorber for zp side
bricks(5).min_x = -1e-3;
bricks(5).min_y = -2e-3;
bricks(5).min_z = 3e-3;
bricks(5).max_x = 11e-3;
bricks(5).max_y = 20e-3;
bricks(5).max_z = 4e-3;
bricks(5).material_type = 5; 

% define absorber for xn side
bricks(6).min_x = -1e-3;
bricks(6).min_y = -2e-3;
bricks(6).min_z = 0;
bricks(6).max_x = 0;
bricks(6).max_y = 20e-3;
bricks(6).max_z = 4e-3;
bricks(6).material_type = 5; 

% define absorber for xp side
bricks(7).min_x = 10e-3;
bricks(7).min_y = -2e-3;
bricks(7).min_z = 0;
bricks(7).max_x = 11e-3;
bricks(7).max_y = 20e-3;
bricks(7).max_z = 4e-3;
bricks(7).material_type = 5; 

% define absorber for yn side
bricks(8).min_x = -1e-3;
bricks(8).min_y = -2e-3;
bricks(8).min_z = 0;
bricks(8).max_x = 11e-3;
bricks(8).max_y = -1e-3;
bricks(8).max_z = 4e-3;
bricks(8).material_type = 5; 

% define absorber for yn side
bricks(9).min_x = -1e-3;
bricks(9).min_y = 19e-3;
bricks(9).min_z = 0;
bricks(9).max_x = 11e-3;
bricks(9).max_y = 20e-3;
bricks(9).max_z = 4e-3;
bricks(9).material_type = 5; 

