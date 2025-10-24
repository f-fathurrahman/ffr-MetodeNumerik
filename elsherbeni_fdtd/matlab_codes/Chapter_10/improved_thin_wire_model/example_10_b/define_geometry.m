disp('defining the problem geometry');

bricks  = [];
spheres = [];
thin_wires = [];

% define a thin wire
thin_wires(1).min_x = 0;
thin_wires(1).min_y = 0;
thin_wires(1).min_z = 0.25e-3;
thin_wires(1).max_x = 0;
thin_wires(1).max_y = 0;
thin_wires(1).max_z = 10e-3;
thin_wires(1).radius = 0.05e-3; 
thin_wires(1).direction = 'z'; 

% define a thin wire
thin_wires(2).min_x = 0;
thin_wires(2).min_y = 0;
thin_wires(2).min_z = -10e-3;
thin_wires(2).max_x = 0;
thin_wires(2).max_y = 0;
thin_wires(2).max_z = -0.25e-3;
thin_wires(2).radius = 0.05e-3; 
thin_wires(2).direction = 'z'; 
