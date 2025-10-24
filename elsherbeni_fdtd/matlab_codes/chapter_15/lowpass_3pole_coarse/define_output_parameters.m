disp('defining output parameters'); 

% figure refresh rate
plotting_step = 100;

% mode of operation
run_simulation = true;
show_material_mesh = true;
show_problem_space = true;

% frequency domain parameters
frequency_domain.start = 0.02e9;
frequency_domain.end   = 10e9;
frequency_domain.step  = 0.02e9;


% a sampled voltage
sampled_voltages(1).min_x = 0;
sampled_voltages(1).min_y = 7.8e-3;
sampled_voltages(1).min_z = 0;
sampled_voltages(1).max_x = 0.4e-3;
sampled_voltages(1).max_y = 9e-3;
sampled_voltages(1).max_z = 0;
sampled_voltages(1).direction = 'xp';
sampled_voltages(1).display_plot = false;

% a sampled voltage
sampled_voltages(2).min_x = 29.6e-3;
sampled_voltages(2).min_y = 7.8e-3;
sampled_voltages(2).min_z = 0;
sampled_voltages(2).max_x = 30e-3;
sampled_voltages(2).max_y = 9e-3;
sampled_voltages(2).max_z = 0;
sampled_voltages(2).direction = 'xn';
sampled_voltages(2).display_plot = false;

% a sampled current
sampled_currents(1).min_x = 0;
sampled_currents(1).min_y = 7.8e-3;
sampled_currents(1).min_z = 0;
sampled_currents(1).max_x = 0.4e-3;
sampled_currents(1).max_y = 9e-3;
sampled_currents(1).max_z = 0;
sampled_currents(1).direction = 'xp';
sampled_currents(1).display_plot = false;

% a sampled current
sampled_currents(2).min_x = 29.6e-3;
sampled_currents(2).min_y = 7.8e-3;
sampled_currents(2).min_z = 0;
sampled_currents(2).max_x = 30e-3;
sampled_currents(2).max_y = 9e-3;
sampled_currents(2).max_z = 0;
sampled_currents(2).direction = 'xn';
sampled_currents(2).display_plot = false;

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

% display problem space parameters
problem_space_display.labels = true;
problem_space_display.axis_at_origin = false;
problem_space_display.axis_outside_domain = true;
problem_space_display.outer_boundaries = true;
problem_space_display.cpml_boundaries = false;
problem_space_display.grid_xn = false;
problem_space_display.grid_yn = false;
problem_space_display.grid_zn = false;
problem_space_display.grid_xp = false;
problem_space_display.grid_yp = false;
problem_space_display.grid_zp = false;
