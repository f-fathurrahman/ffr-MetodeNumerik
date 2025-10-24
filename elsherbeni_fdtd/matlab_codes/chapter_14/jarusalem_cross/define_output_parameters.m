disp('defining output parameters'); 

% figure refresh rate
plotting_step = 20;

% mode of operation
run_simulation = true;
show_material_mesh = true;
show_problem_space = true;

% frequency domain parameters
frequency_domain.start = 0.2;
frequency_domain.end   = 4;
frequency_domain.step  = 0.2;

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
