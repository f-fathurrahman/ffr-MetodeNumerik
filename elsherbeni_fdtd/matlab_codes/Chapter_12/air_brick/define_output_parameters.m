disp('defining output parameters'); 

% figure refresh rate
plotting_step = 5;

% mode of operation
run_simulation = true;
show_material_mesh = false;
show_problem_space = false;

% frequency domain parameters
frequency_domain.start = 0.2;
frequency_domain.end   = 4;
frequency_domain.step  = 0.2;

% define animation
% field_type shall be 'e' or 'h'
% plane cut shall be 'xy', 'yz', or 'zx'
% component shall be 'x', 'y', 'z', or 'm'
% an animation
animation(1).field_type = 'e';
animation(1).component = 'm';
animation(1).plane_cut(1).type = 'yz';
animation(1).plane_cut(1).position  = 10e-3;
animation(1).display_grid = false;
animation(1).display_objects = true;
animation(1).save_movie = true;
animation(1).view_angles = [55 15];
animation(1).zoom = 0.8;
animation(1).enable = true;

% display problem space parameters
problem_space_display.labels = true;
problem_space_display.axis_at_origin = false;
problem_space_display.axis_outside_domain = true;
problem_space_display.outer_boundaries = true;
problem_space_display.cpml_boundaries = true;
problem_space_display.grid_xn = false;
problem_space_display.grid_yn = false;
problem_space_display.grid_zn = true;
problem_space_display.grid_xp = true;
problem_space_display.grid_yp = true;
problem_space_display.grid_zp = false;
