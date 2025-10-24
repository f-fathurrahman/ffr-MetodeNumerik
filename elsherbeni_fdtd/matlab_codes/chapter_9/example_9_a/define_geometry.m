disp('defining the problem geometry');

bricks  = [];
spheres = [];

% define substrate
bricks(1).min_x = -0.787e-3;
bricks(1).min_y = 0;
bricks(1).min_z = 0;
bricks(1).max_x = 0;
bricks(1).max_y = 40e-3;
bricks(1).max_z = 40e-3;
bricks(1).material_type = 4; 

bricks(2).min_x = 0;
bricks(2).min_y = 0;
bricks(2).min_z = 24e-3;
bricks(2).max_x = 0;
bricks(2).max_y = 28.4e-3;
bricks(2).max_z = 26.4e-3;
bricks(2).material_type = 2; 

bricks(3).min_x = 0;
bricks(3).min_y = 16e-3;
bricks(3).min_z = 30e-3;
bricks(3).max_x = 0;
bricks(3).max_y = 28.4e-3;
bricks(3).max_z = 32.4e-3;
bricks(3).material_type = 2; 

bricks(4).min_x = 0;
bricks(4).min_y = 26e-3;
bricks(4).min_z = 8.4e-3;
bricks(4).max_x = 0;
bricks(4).max_y = 28.4e-3;
bricks(4).max_z = 32.4e-3;
bricks(4).material_type = 2; 

bricks(5).min_x = 0;
bricks(5).min_y = 20.8e-3;
bricks(5).min_z = 16e-3;
bricks(5).max_x = 0;
bricks(5).max_y = 23.2e-3;
bricks(5).max_z = 32.4e-3;
bricks(5).material_type = 2; 

bricks(6).min_x = -0.787e-3;
bricks(6).min_y = 16e-3;
bricks(6).min_z = 30e-3;
bricks(6).max_x = 0;
bricks(6).max_y = 16e-3;
bricks(6).max_z = 32.4e-3;
bricks(6).material_type = 2; 

bricks(7).min_x = -0.787e-3;
bricks(7).min_y = 0;
bricks(7).min_z = 0;
bricks(7).max_x = -0.787e-3;
bricks(7).max_y = 16e-3;
bricks(7).max_z = 40e-3;
bricks(7).material_type = 2; 
