disp('defining output parameters'); 

% figure refresh rate
plotting_step = 20;

% mode of operation
run_simulation = true;
show_material_mesh = true;
show_problem_space = true;

% frequency domain parameters
frequency_domain.start = 0.2;
frequency_domain.end   = 20;
frequency_domain.step  = 0.2;

sampled_electric_fields(1).component = 'm';
sampled_electric_fields(1).x = 5e-3;
sampled_electric_fields(1).y = 5e-3;
sampled_electric_fields(1).z = 1e-2;
sampled_electric_fields(1).display_plot = false;


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
