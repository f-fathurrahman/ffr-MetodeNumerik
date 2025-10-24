disp('defining output parameters'); 

% figure refresh rate
plotting_step = 100;

% mode of operation
run_simulation = true;
show_material_mesh = true;
show_problem_space = true;

% frequency domain parameters
frequency_domain.start = 0.2e9;
frequency_domain.end   = 12e9;
frequency_domain.step  = 0.02e9;

% a sampled voltage
sampled_voltages(1).min_x = 29e-3;
sampled_voltages(1).min_y = 10e-3;
sampled_voltages(1).min_z = 0;
sampled_voltages(1).max_x = 31e-3;
sampled_voltages(1).max_y = 10e-3;
sampled_voltages(1).max_z = 1.9e-3;
sampled_voltages(1).direction = 'zp';
sampled_voltages(1).display_plot = false;

% a sampled current
sampled_currents(1).min_x = 29e-3;
sampled_currents(1).min_y = 10e-3;
sampled_currents(1).min_z = 0;
sampled_currents(1).max_x = 31e-3;
sampled_currents(1).max_y = 10e-3;
sampled_currents(1).max_z = 1.9e-3;
sampled_currents(1).direction = 'zp';
sampled_currents(1).display_plot = false;

% a port
ports(1).sampled_voltage_index = 1;
ports(1).sampled_current_index = 1;
ports(1).impedance = 50;
ports(1).is_source_port = true;
 
% % display problem space parameters
% problem_space_display.labels = true;
% problem_space_display.axis_at_origin = false;
% problem_space_display.axis_outside_domain = true;
% problem_space_display.outer_boundaries = true;
% problem_space_display.cpml_boundaries = true;
% problem_space_display.grid_xn = false;
% problem_space_display.grid_yn = false;
% problem_space_display.grid_zn = true;
% problem_space_display.grid_xp = true;
% problem_space_display.grid_yp = true;
% problem_space_display.grid_zp = false;
