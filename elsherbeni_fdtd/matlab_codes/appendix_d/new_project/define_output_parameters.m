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


% a sampled voltage
sampled_voltages(1).min_x = 0.1;
sampled_voltages(1).min_y = 0.1;
sampled_voltages(1).min_z = 0;
sampled_voltages(1).max_x = 0.1;
sampled_voltages(1).max_y = 0.1;
sampled_voltages(1).max_z = 0.1;
sampled_voltages(1).direction = 'zp';
sampled_voltages(1).display_plot = true;

% a sampled voltage
sampled_voltages(2).min_x = 0.2;
sampled_voltages(2).min_y = 0.2;
sampled_voltages(2).min_z = 0;
sampled_voltages(2).max_x = 0.2;
sampled_voltages(2).max_y = 0.2;
sampled_voltages(2).max_z = 0.2;
sampled_voltages(2).direction = 'zp';
sampled_voltages(2).display_plot = true;


% a sampled current
sampled_currents(1).min_x = 0.1;
sampled_currents(1).min_y = 0.1;
sampled_currents(1).min_z = 0;
sampled_currents(1).max_x = 0.1;
sampled_currents(1).max_y = 0.1;
sampled_currents(1).max_z = 0.1;
sampled_currents(1).direction = 'zp';
sampled_currents(1).display_plot = true;

% a sampled current
sampled_currents(2).min_x = 0.2;
sampled_currents(2).min_y = 0.2;
sampled_currents(2).min_z = 0;
sampled_currents(2).max_x = 0.2;
sampled_currents(2).max_y = 0.2;
sampled_currents(2).max_z = 0.2;
sampled_currents(2).direction = 'zp';
sampled_currents(2).display_plot = true;


% a port
ports(1).sampled_voltage_index = 1;
ports(1).sampled_current_index = 1;
ports(1).impedance = 50;
ports(1).is_source_port = true;

% a port
ports(2).sampled_voltage_index = 2;
ports(2).sampled_current_index = 2;
ports(2).impedance = 50;
ports(2).is_source_port = false;


% define animation
% field_type shall be 'e' or 'h'
% plane cut shall be 'xy', 'yz', or 'zx'
% component shall be 'x', 'y', 'z', or 'm'
% an animation
animation(1).field_type = 'e';
animation(1).component = 'm';
animation(1).plane_cut(1).type = 'xy';
animation(1).plane_cut(1).position  = 0;
animation(1).display_grid = true;
animation(1).display_objects = true;
animation(1).save_movie = true;
animation(1).view_angles = [40 30];
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
