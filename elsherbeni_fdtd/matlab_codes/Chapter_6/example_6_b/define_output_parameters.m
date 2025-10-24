disp('defining output parameters');

sampled_electric_fields = [];
sampled_magnetic_fields = [];
sampled_voltages = [];
sampled_currents = [];

% figure refresh rate
plotting_step = 100;

% mode of operation
run_simulation = true;
show_material_mesh = true;
show_problem_space = true;

% frequency domain parameters
frequency_domain.start = 2e7;
frequency_domain.end   = 8e9;
frequency_domain.step  = 2e7;

% define sampled voltages
sampled_voltages(1).min_x = 4e-3;
sampled_voltages(1).min_y = 0;
sampled_voltages(1).min_z = 0;
sampled_voltages(1).max_x = 5.8e-3;
sampled_voltages(1).max_y = 0.4e-3;
sampled_voltages(1).max_z = 1e-3;
sampled_voltages(1).direction = 'zp';
sampled_voltages(1).display_plot = false;

% define sampled voltages
sampled_voltages(2).min_x = 4.8e-3;
sampled_voltages(2).min_y = 17.6e-3;
sampled_voltages(2).min_z = 0;
sampled_voltages(2).max_x = 5.2e-3;
sampled_voltages(2).max_y = 18e-3;
sampled_voltages(2).max_z = 1e-3;
sampled_voltages(2).direction = 'zp';
sampled_voltages(2).display_plot = false;

% define sampled currents
sampled_currents(1).min_x = 4e-3;
sampled_currents(1).min_y = 0;
sampled_currents(1).min_z = 0.4e-3;
sampled_currents(1).max_x = 5.8e-3;
sampled_currents(1).max_y = 0.4e-3;
sampled_currents(1).max_z = 0.6e-3;
sampled_currents(1).direction = 'zp';
sampled_currents(1).display_plot = false;

% define sampled currents
sampled_currents(2).min_x = 4.8e-3;
sampled_currents(2).min_y = 17.6e-3;
sampled_currents(2).min_z = 0.4e-3;
sampled_currents(2).max_x = 5.2e-3;
sampled_currents(2).max_y = 18e-3;
sampled_currents(2).max_z = 0.6e-3;
sampled_currents(2).direction = 'zp';
sampled_currents(2).display_plot = false;

% define ports
ports(1).sampled_voltage_index = 1;
ports(1).sampled_current_index = 1;
ports(1).impedance = 50;
ports(1).is_source_port = true;

ports(2).sampled_voltage_index = 2;
ports(2).sampled_current_index = 2;
ports(2).impedance = 100;
ports(2).is_source_port = false;


% define animation
% field_type shall be 'e' or 'h'
% plane cut shall be 'xy', yz, or zx
% component shall be 'x', 'y', 'z', or 'm;
animation(1).field_type = 'e';
animation(1).component = 'x';
animation(1).plane_cut(1).type = 'xy';
animation(1).plane_cut(1).position  = 1e-3;
animation(1).enable = false;
animation(1).display_grid = false;
animation(1).display_objects = true;

% display problem space parameters
problem_space_display.labels = false;
problem_space_display.axis_at_origin = false;
problem_space_display.axis_outside_domain = true;
problem_space_display.grid_xn = false;
problem_space_display.grid_xp = false;
problem_space_display.grid_yn = false;
problem_space_display.grid_yp = false;
problem_space_display.grid_zn = false;
problem_space_display.grid_zp = false;
problem_space_display.outer_boundaries = false;
problem_space_display.cpml_boundaries = false;

