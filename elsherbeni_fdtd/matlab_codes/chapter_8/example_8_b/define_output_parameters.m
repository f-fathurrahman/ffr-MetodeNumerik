disp('defining output parameters');

sampled_electric_fields = [];
sampled_magnetic_fields = [];
sampled_voltages = [];
sampled_currents = [];
ports = [];

% figure refresh rate
plotting_step = 5;

% mode of operation
run_simulation = true;
show_material_mesh = true;
show_problem_space = true;

% frequency domain parameters
frequency_domain.start = 20e6;
frequency_domain.end   = 10e9;
frequency_domain.step  = 20e6;

% define sampled voltages
sampled_voltages(1).min_x = -15*dx;
sampled_voltages(1).min_y = -25*dy;
sampled_voltages(1).min_z = 0;
sampled_voltages(1).max_x = -9*dx;
sampled_voltages(1).max_y = -25*dy;
sampled_voltages(1).max_z = 3*dz;
sampled_voltages(1).direction = 'zp';
sampled_voltages(1).display_plot = false;

% define sampled voltages
sampled_voltages(2).min_x = 9*dx;
sampled_voltages(2).min_y = -25*dy;
sampled_voltages(2).min_z = 0;
sampled_voltages(2).max_x = 15*dx;
sampled_voltages(2).max_y = -25*dy;
sampled_voltages(2).max_z = 3*dz;
sampled_voltages(2).direction = 'zp';
sampled_voltages(2).display_plot = false;

% define sampled voltages
sampled_voltages(3).min_x = 9*dx;
sampled_voltages(3).min_y = 25*dy;
sampled_voltages(3).min_z = 0;
sampled_voltages(3).max_x = 15*dx;
sampled_voltages(3).max_y = 25*dy;
sampled_voltages(3).max_z = 3*dz;
sampled_voltages(3).direction = 'zp';
sampled_voltages(3).display_plot = false;

% define sampled voltages
sampled_voltages(4).min_x = -15*dx;
sampled_voltages(4).min_y = 25*dy;
sampled_voltages(4).min_z = 0;
sampled_voltages(4).max_x = -9*dx;
sampled_voltages(4).max_y = 25*dy;
sampled_voltages(4).max_z = 3*dz;
sampled_voltages(4).direction = 'zp';
sampled_voltages(4).display_plot = false;

% define sampled currents
sampled_currents(1).min_x = -15*dx;
sampled_currents(1).min_y = -25*dy;
sampled_currents(1).min_z = 2*dz;
sampled_currents(1).max_x = -9*dx;
sampled_currents(1).max_y = -25*dy;
sampled_currents(1).max_z = 2*dz;
sampled_currents(1).direction = 'zp';
sampled_currents(1).display_plot = false;

% define sampled currents
sampled_currents(2).min_x = 9*dx;
sampled_currents(2).min_y = -25*dy;
sampled_currents(2).min_z = 2*dz;
sampled_currents(2).max_x = 15*dx;
sampled_currents(2).max_y = -25*dy;
sampled_currents(2).max_z = 2*dz;
sampled_currents(2).direction = 'zp';
sampled_currents(2).display_plot = false;

% define sampled currents
sampled_currents(3).min_x = 9*dx;
sampled_currents(3).min_y = 25*dy;
sampled_currents(3).min_z = 2*dz;
sampled_currents(3).max_x = 15*dx;
sampled_currents(3).max_y = 25*dy;
sampled_currents(3).max_z = 2*dz;
sampled_currents(3).direction = 'zp';
sampled_currents(3).display_plot = false;

% define sampled currents
sampled_currents(4).min_x = -15*dx;
sampled_currents(4).min_y = 25*dy;
sampled_currents(4).min_z = 2*dz;
sampled_currents(4).max_x = -9*dx;
sampled_currents(4).max_y = 25*dy;
sampled_currents(4).max_z = 2*dz;
sampled_currents(4).direction = 'zp';
sampled_currents(4).display_plot = false;

% define ports
ports(1).sampled_voltage_index = 1;
ports(1).sampled_current_index = 1;
ports(1).impedance = 50;
ports(1).is_source_port = true;

ports(2).sampled_voltage_index = 2;
ports(2).sampled_current_index = 2;
ports(2).impedance = 50;
ports(2).is_source_port = false;

ports(3).sampled_voltage_index = 3;
ports(3).sampled_current_index = 3;
ports(3).impedance = 50;
ports(3).is_source_port = false;

ports(4).sampled_voltage_index = 4;
ports(4).sampled_current_index = 4;
ports(4).impedance = 50;
ports(4).is_source_port = false;

% define animation
% field_type shall be 'e' or 'h'
% plane cut shall be 'xy', yz, or zx
% component shall be 'x', 'y', 'z', or 'm;
animation(1).field_type = 'e';
animation(1).component = 'm';
animation(1).plane_cut(1).type = 'xy';
animation(1).plane_cut(1).position  = 2*dz;
animation(1).enable = true;
animation(1).display_grid = false;
animation(1).display_objects = true;

% display problem space parameters
problem_space_display.labels = false;
problem_space_display.axis_at_origin = false;
problem_space_display.axis_outside_domain = false;
problem_space_display.grid_xn = false;
problem_space_display.grid_xp = false;
problem_space_display.grid_yn = false;
problem_space_display.grid_yp = false;
problem_space_display.grid_zn = false;
problem_space_display.grid_zp = false;
problem_space_display.outer_boundaries = false;
problem_space_display.cpml_boundaries = false;
