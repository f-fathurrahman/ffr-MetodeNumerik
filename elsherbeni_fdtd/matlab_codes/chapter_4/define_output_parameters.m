disp('defining output parameters');

sampled_electric_fields = [];
sampled_magnetic_fields = [];
sampled_voltages = [];
sampled_currents = [];

% figure refresh rate
plotting_step = 10;

% mode of operation
run_simulation = true;
show_material_mesh = true;
show_problem_space = true;

% define sampled electric fields
% component: vector component 'x','y','z', or magnitude 'm'
% display_plot = true, in order to plot field during simulation 
sampled_electric_fields(1).x = 30*dx;
sampled_electric_fields(1).y = 30*dy;
sampled_electric_fields(1).z = 10*dz;
sampled_electric_fields(1).component = 'x';
sampled_electric_fields(1).display_plot = true;

% define sampled magnetic fields
% component: vector component 'x','y','z', or magnitude 'm'
% display_plot = true, in order to plot field during simulation 
sampled_magnetic_fields(1).x = 30*dx;
sampled_magnetic_fields(1).y = 30*dy;
sampled_magnetic_fields(1).z = 10*dz;
sampled_magnetic_fields(1).component = 'm';
sampled_magnetic_fields(1).display_plot = true;

% define sampled voltages
sampled_voltages(1).min_x = 5.0e-3;
sampled_voltages(1).min_y = 0;
sampled_voltages(1).min_z = 0;
sampled_voltages(1).max_x = 5.0e-3;
sampled_voltages(1).max_y = 2.0e-3;
sampled_voltages(1).max_z = 4.0e-3;
sampled_voltages(1).direction = 'zp';
sampled_voltages(1).display_plot = true;

% define sampled currents
sampled_currents(1).min_x = 5.0e-3;
sampled_currents(1).min_y = 0;
sampled_currents(1).min_z = 4.0e-3;
sampled_currents(1).max_x = 5.0e-3;
sampled_currents(1).max_y = 2.0e-3;
sampled_currents(1).max_z = 4.0e-3;
sampled_currents(1).direction = 'xp';
sampled_currents(1).display_plot = true;

% display problem space parameters
problem_space_display.labels = true;
problem_space_display.axis_at_origin = false;
problem_space_display.axis_outside_domain = true;
problem_space_display.grid_xn = false;
problem_space_display.grid_xp = true;
problem_space_display.grid_yn = false;
problem_space_display.grid_yp = true;
problem_space_display.grid_zn = true;
problem_space_display.grid_zp = false;
problem_space_display.outer_boundaries = true;
problem_space_display.cpml_boundaries = true;

% define animation
% field_type shall be 'e' or 'h'
% plane cut shall be 'xy', yz, or zx
% component shall be 'x', 'y', 'z', or 'm;
animation(1).field_type = 'e';
animation(1).component = 'm';
animation(1).plane_cut(1).type = 'xy';
animation(1).plane_cut(1).position  = 0;
animation(1).plane_cut(2).type = 'yz';
animation(1).plane_cut(2).position  = 0;
animation(1).plane_cut(3).type = 'zx';
animation(1).plane_cut(3).position  = 0;
animation(1).enable = true;
animation(1).display_grid = false;
animation(1).display_objects = true;

animation(2).field_type = 'h';
animation(2).component = 'x';
animation(2).plane_cut(1).type = 'xy';
animation(2).plane_cut(1).position  = -5;
animation(2).plane_cut(2).type = 'xy';
animation(2).plane_cut(2).position  = 5;
animation(2).enable = true;
animation(2).display_grid = true;
animation(2).display_objects = true;
